package wdl4s.cwl

import shapeless.syntax.singleton._
import shapeless._
import cats.syntax.foldable._
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
import wdl4s.wom.callable.Callable.{OutputDefinition, RequiredInputDefinition}
import wdl4s.wom.callable.{Callable, TaskDefinition}
import wdl4s.wom.executable.Executable
import wdl4s.wom.expression.{WomExpression, PlaceholderWomExpression}
import wdl4s.wom.graph._

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

                            //TODO: Fix this when CwlAny parses correctly
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
                            permanentFailCodes: Option[Array[Int]] = None) {

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

  def graphNodes: ErrorOr[Set[GraphNode]] = {
    CallNode.
      callWithInputs(id.getOrElse("this is a made up call node name"), taskDefinition, Map.empty, Set.empty).
      map(cwi => Set.empty[GraphNode] ++ cwi.inputs + cwi.call)
  }

  def asCwl = Coproduct[Cwl](this)
}

object CommandLineTool {

  type StringOrExpression = ECMAScriptExpression :+: String :+: CNil

  type BaseCommand = String :+: Array[String] :+: CNil

  type Argument = ECMAScriptExpression :+: CommandLineBinding :+: String :+: CNil
}

