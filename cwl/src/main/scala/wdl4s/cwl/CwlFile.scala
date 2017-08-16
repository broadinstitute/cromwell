package wdl4s.cwl

import shapeless.syntax.singleton._
import cats.syntax.foldable._
import cats.syntax.traverse._
import cats.instances.list._
import cats.syntax.option._
import shapeless.{:+:, CNil, Poly1, Witness}
import CwlType._
import shapeless.syntax.singleton._
import CwlVersion._
import cats.data.Validated._
import lenthall.validation.ErrorOr._
import wdl4s.cwl.CommandLineTool.{BaseCommand, StringOrExpression}
import wdl4s.cwl.CwlType.CwlType
import wdl4s.wdl.{RuntimeAttributes, WdlExpression}
import wdl4s.wdl.command.CommandPart
import wdl4s.wdl.types._
import wdl4s.wom.callable.Callable.{OutputDefinition, RequiredInputDefinition}
import wdl4s.wom.callable.{Callable, TaskDefinition, WorkflowDefinition}
import wdl4s.wom.executable.Executable
import wdl4s.wom.expression.{WomExpression, PlaceholderWomExpression}
import wdl4s.wom.graph.GraphNodePort.{GraphNodeOutputPort, OutputPort}
import wdl4s.wom.graph._

sealed trait CwlFile {

  val cwlVersion: Option[CwlVersion]
}

case class Workflow(
                     cwlVersion: Option[CwlVersion] = Option(CwlVersion.Version1),
                     `class`: Workflow.`class`.type = Workflow.`class`,
                     inputs: Array[InputParameter] = Array.empty,
                     outputs: Array[WorkflowOutputParameter] = Array.empty,
                     steps: Array[WorkflowStep]) extends CwlFile {

  def womExecutable(cwlFileMap: Map[String, CwlFile]): ErrorOr[Executable] = womDefinition(cwlFileMap) map Executable.apply

  def outputsTypeMap(cwlFileMap: Map[String, CwlFile]): WdlTypeMap = steps.foldLeft(Map.empty[String, WdlType]) {
    // Not implemented as a `foldMap` because there is no semigroup instance for `WdlType`s.  `foldMap` doesn't know that
    // we don't need a semigroup instance since the map keys should be unique and therefore map values would never need
    // to be combined under the same key.
    (acc, s) => acc ++ s.typedOutputs(cwlFileMap)
  }

  lazy val stepById: Map[String, WorkflowStep] = steps.map(ws => ws.id -> ws).toMap

  def womGraph(cwlFileMap: Map[String, CwlFile]): ErrorOr[Graph] = {

    def cwlTypeForInputParameter(input: InputParameter): Option[CwlType] = input.`type`.flatMap(_.select[CwlType])

    def wdlTypeForInputParameter(input: InputParameter): Option[WdlType] = {
      cwlTypeForInputParameter(input) map cwlTypeToWdlType
    }

    val typeMap: WdlTypeMap =
      outputsTypeMap(cwlFileMap) ++
        // Note this is only looking at the workflow inputs and not recursing into steps, because our current thinking
        // is that in CWL graph inputs can only be defined at the workflow level.  It's possible that's not actually
        // correct, but that's the assumption being made here.
        inputs.toList.flatMap { i =>
          wdlTypeForInputParameter(i).map(i.id -> _).toList
        }.toMap

    val graphFromInputs: Set[RequiredGraphInputNode] = inputs.map {
      input => RequiredGraphInputNode(input.id, wdlTypeForInputParameter(input).get)
    }.toSet

    val workflowInputs: Map[String, GraphNodeOutputPort] =
      graphFromInputs.map {
        workflowInput =>
          workflowInput.name -> GraphNodeOutputPort(workflowInput.name, workflowInput.womType, workflowInput)
      }.toMap

    val graphFromSteps =
      steps.
        toList.
        foldLeft(Set.empty[GraphNode])(
          (nodes, step) => step.callWithInputs(typeMap, cwlFileMap, this, nodes, workflowInputs))

    val graphFromOutputs: ErrorOr[Set[GraphNode]] =
      outputs.toList.traverse[ErrorOr, GraphNode] {
        output =>

          val wdlType = cwlTypeToWdlType(output.`type`.flatMap(_.select[CwlType]).get)

          def lookupOutputSource(source: String): ErrorOr[OutputPort] =
            graphFromSteps.
              flatMap(_.outputPorts).
              find(_.name == source).
              toValidNel(s"unable to find upstream port corresponding to $source")

          lookupOutputSource(output.outputSource.flatMap(_.select[String]).get).
            map(PortBasedGraphOutputNode(output.id, wdlType, _))
      }.map(_.toSet)

    graphFromOutputs.flatMap(outputs => Graph.validateAndConstruct(graphFromSteps ++ graphFromInputs ++ outputs))
  }

  def womDefinition(cwlFileMap: Map[String, CwlFile]): ErrorOr[WorkflowDefinition] = {
    val name: String = "workflow Id"
    val meta: Map[String, String] = Map.empty
    val paramMeta: Map[String, String] = Map.empty
    val declarations: List[(String, WomExpression)] = List.empty

    womGraph(cwlFileMap).map(graph =>
      WorkflowDefinition(
        name,
        graph,
        meta,
        paramMeta,
        declarations
      )
    )
  }
}


object Workflow {
  val `class`: Witness.`"Workflow"`.T = "Workflow".narrow
}

/**
  *
  * @param inputs
  * @param outputs
  * @param `class` This _should_ always be "CommandLineTool," however the spec does not -er- specify this.
  * @param id
  * @param requirements
  * @param hints
  * @param label
  * @param doc
  * @param cwlVersion
  * @param baseCommand
  * @param arguments
  * @param stdin
  * @param stderr
  * @param stdout
  * @param successCodes
  * @param temporaryFailCodes
  * @param permanentFailCodes
  */
case class CommandLineTool(
                            inputs: Array[CommandInputParameter] = Array.empty,
                            outputs: Array[CommandOutputParameter] = Array.empty,
                            `class`: Witness.`"CommandLineTool"`.T = "CommandLineTool".narrow,
                            id: Option[String] = None,
                            requirements: Option[Array[Requirement]] = None,
                            //hints: Option[Array[CwlAny]] = None,
                            hints: Option[Array[Map[String, String]]] = None,
                            label: Option[String] = None,
                            doc: Option[String] = None,
                            cwlVersion: Option[CwlVersion] = Option(CwlVersion.Version1),
                            baseCommand: Option[BaseCommand] = None,
                            arguments: Option[Array[CommandLineTool.Argument]] = None,
                            stdin: Option[StringOrExpression] = None,
                            stderr: Option[StringOrExpression] = None,
                            stdout: Option[StringOrExpression] = None,
                            successCodes: Option[Array[Int]] = None,
                            temporaryFailCodes: Option[Array[Int]] = None,
                            permanentFailCodes: Option[Array[Int]] = None) extends CwlFile {

  def womExecutable: ErrorOr[Executable] =
    Valid(Executable(taskDefinition))


  object BaseCommandToString extends Poly1 {
    implicit def one = at[String] {
      identity
    }

    implicit def many = at[Array[String]] {
      _.mkString(" && ")
    }
  }

  object ArgumentToId extends Poly1 {
    implicit def ecmaScript = at[ECMAScriptExpression] {
      _.value
    }

    implicit def commandLineBinding = at[CommandLineBinding] { _ => "" }

    implicit def string = at[String] {
      identity
    }
  }

  /**
    * This is used in place of the id when id is None.
    *
    * @return
    */
  def taskDefinitionId: String =
    baseCommand.map(_.fold(BaseCommandToString)).getOrElse(
      arguments.map(_.map(_.fold(ArgumentToId)).mkString(" ")).get)

  def taskDefinition: TaskDefinition = {

    val id = this.id.getOrElse(taskDefinitionId)

    val commandTemplate: Seq[CommandPart] = baseCommand.get.fold(BaseCommandToCommandParts)

    val runtimeAttributes: RuntimeAttributes = RuntimeAttributes(Map.empty[String, WdlExpression])

    val meta: Map[String, String] = Map.empty
    val parameterMeta: Map[String, String] = Map.empty

    //TODO: This output does _not_ capture expressions from the output.outputBinding
    //The implementation must include the expression evaluation pieces as detailed in:
    //http://www.commonwl.org/v1.0/CommandLineTool.html#CommandOutputBinding
    val outputs: Set[Callable.OutputDefinition] = this.outputs.map {
      output =>
        val wdlType = output.`type`.flatMap(_.select[CwlType]).map(cwlTypeToWdlType).get //<-- here be `get` dragons
        OutputDefinition(output.id, wdlType, PlaceholderWomExpression(Set.empty, wdlType))
    }.toSet

    val inputs: Set[_ <: Callable.InputDefinition] =
      this.inputs.map { cip =>
        val tpe = cip.`type`.flatMap(_.select[CwlType]).map(cwlTypeToWdlType).get

        //TODO: This id includes the filename, which makes assigning input values more laborious
        //We should consider dropping filenames for _all_ ids, as long as we can guarantee uniqueness
        val inputId = cip.id
        RequiredInputDefinition(inputId, tpe)
      }.toSet

    val declarations: List[(String, WomExpression)] = List.empty

    TaskDefinition(
      id,
      commandTemplate,
      runtimeAttributes,
      meta,
      parameterMeta,
      outputs,
      inputs,
      declarations
    )
  }

  def graphNodes: Set[GraphNode] = {

    //need to gather up the step outputs and pass them into Call with inputs


    val cwi = CallNode.callWithInputs(id.getOrElse("this is a made up call node name"), taskDefinition, Map.empty)

    Set.empty[GraphNode] ++ cwi.inputs + cwi.call
  }
}

object CommandLineTool {

  type StringOrExpression = ECMAScriptExpression :+: String :+: CNil

  type BaseCommand = String :+: Array[String] :+: CNil

  type Argument = ECMAScriptExpression :+: CommandLineBinding :+: String :+: CNil
}

