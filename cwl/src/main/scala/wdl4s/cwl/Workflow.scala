package wdl4s.cwl

import cats.data.Validated._
import cats.instances.list._
import cats.syntax.option._
import cats.syntax.traverse._
import CwlType._
import shapeless._
import shapeless.syntax.singleton._
import CwlVersion._
import cats.data.Validated._
import lenthall.validation.ErrorOr._
import shapeless._
import wdl4s.cwl.CwlType.CwlType
import wdl4s.cwl.CwlVersion._
import wdl4s.wdl.types._
import wdl4s.wom.callable.WorkflowDefinition
import wdl4s.wom.executable.Executable
import wdl4s.wom.expression.WomExpression
import wdl4s.wom.graph.GraphNodePort.{GraphNodeOutputPort, OutputPort}
import wdl4s.wom.graph._

case class Workflow private(
                     cwlVersion: Option[CwlVersion],
                     `class`: Witness.`"Workflow"`.T,
                     inputs: Array[InputParameter],
                     outputs: Array[WorkflowOutputParameter],
                     steps: Array[WorkflowStep]) {

  def womExecutable: ErrorOr[Executable] = womDefinition map Executable.apply

  val fileNames: List[String] = steps.toList.flatMap(_.run.select[String].toList)

  def outputsTypeMap: WdlTypeMap = steps.foldLeft(Map.empty[String, WdlType]) {
    // Not implemented as a `foldMap` because there is no semigroup instance for `WdlType`s.  `foldMap` doesn't know that
    // we don't need a semigroup instance since the map keys should be unique and therefore map values would never need
    // to be combined under the same key.
    (acc, s) => acc ++ s.typedOutputs
  }

  lazy val stepById: Map[String, WorkflowStep] = steps.map(ws => ws.id -> ws).toMap

  def womGraph: ErrorOr[Graph] = {

    def cwlTypeForInputParameter(input: InputParameter): Option[CwlType] = input.`type`.flatMap(_.select[CwlType])

    def wdlTypeForInputParameter(input: InputParameter): Option[WdlType] = {
      cwlTypeForInputParameter(input) map cwlTypeToWdlType
    }

    val typeMap: WdlTypeMap =
      outputsTypeMap ++
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
          workflowInput.name -> workflowInput.singleOutputPort
      }.toMap

    val graphFromSteps: Set[GraphNode] =
      steps.
        toList.
        foldLeft(Set.empty[GraphNode] ++ graphFromInputs)(
          (nodes, step) => step.callWithInputs(typeMap,  this, nodes, workflowInputs))

    val graphFromOutputs: ErrorOr[Set[GraphNode]] =
      outputs.toList.traverse[ErrorOr, GraphNode] {
        output =>

          val wdlType = cwlTypeToWdlType(output.`type`.flatMap(_.select[CwlType]).get)

          def lookupOutputSource(outputId: WorkflowOutputId): ErrorOr[OutputPort] = {
            (for {
              call <- graphFromSteps.collectFirst { case callNode: CallNode if callNode.name == outputId.stepId => callNode }
              output <- call.outputPorts.find(_.name == outputId.outputId)
            } yield output).toValidNel(s"unable to find upstream port corresponding to ${outputId.stepId}/${outputId.outputId}")
          }

          lookupOutputSource(WorkflowOutputId(output.outputSource.flatMap(_.select[String]).get)).
            map(PortBasedGraphOutputNode(output.id, wdlType, _))
      }.map(_.toSet)

    graphFromOutputs.flatMap(outputs => Graph.validateAndConstruct(graphFromSteps ++ graphFromInputs ++ outputs))
  }

  def womDefinition: ErrorOr[WorkflowDefinition] = {
    // TODO: need to find a way to get a meaningful name here
    val name: String = "MyCwlWorkflow"
    val meta: Map[String, String] = Map.empty
    val paramMeta: Map[String, String] = Map.empty
    val declarations: List[(String, WomExpression)] = List.empty

    womGraph.map(graph =>
      WorkflowDefinition(
        name,
        graph,
        meta,
        paramMeta,
        declarations
      )
    )
  }

  def asCwl = Coproduct[Cwl](this)
}
object Workflow {

  def apply(cwlVersion: Option[CwlVersion] = Option(CwlVersion.Version1),
            inputs: Array[InputParameter] = Array.empty,
            outputs: Array[WorkflowOutputParameter] = Array.empty,
            steps: Array[WorkflowStep] = Array.empty): Workflow  =
              Workflow(cwlVersion, "Workflow".narrow, inputs, outputs, steps)
}
