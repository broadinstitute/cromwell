package cwl

import cats.data.NonEmptyList
import cats.instances.list._
import cats.syntax.either._
import cats.syntax.traverse._
import lenthall.Checked
import lenthall.validation.ErrorOr._
import shapeless._
import shapeless.syntax.singleton._
import CwlType.CwlType
import CwlVersion._
import wom.callable.WorkflowDefinition
import wom.executable.Executable
import wom.expression.{PlaceholderWomExpression, WomExpression}
import wom.graph.GraphNodePort.{GraphNodeOutputPort, OutputPort}
import wom.graph._
import wom.types.WdlType

case class Workflow private(
                     cwlVersion: Option[CwlVersion],
                     `class`: Witness.`"Workflow"`.T,
                     inputs: Array[InputParameter],
                     outputs: Array[WorkflowOutputParameter],
                     steps: Array[WorkflowStep]) {
  def womExecutable(inputFile: Option[String] = None): Checked[Executable] = {
    CwlExecutableValidation.builWomExecutable(womDefinition, inputFile)
  }

  val fileNames: List[String] = steps.toList.flatMap(_.run.select[String].toList)

  def outputsTypeMap: WdlTypeMap = steps.foldLeft(Map.empty[String, WdlType]) {
    // Not implemented as a `foldMap` because there is no semigroup instance for `WdlType`s.  `foldMap` doesn't know that
    // we don't need a semigroup instance since the map keys should be unique and therefore map values would never need
    // to be combined under the same key.
    (acc, s) => acc ++ s.typedOutputs
  }

  lazy val stepById: Map[String, WorkflowStep] = steps.map(ws => ws.id -> ws).toMap

  def womGraph: Checked[Graph] = {

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

    val graphFromInputs: Set[ExternalGraphInputNode] = inputs.map {
      // TODO WOM: need to be able to transform this default value to a WomExpression
      case inputParameter if inputParameter.default.isDefined =>
        val parsedInputId = WorkflowInputId(inputParameter.id).inputId
        val womType = wdlTypeForInputParameter(inputParameter).get

        OptionalGraphInputNodeWithDefault(WomIdentifier(parsedInputId), womType, PlaceholderWomExpression(Set.empty, womType))
      case input =>
        val parsedInputId = WorkflowInputId(input.id).inputId

        RequiredGraphInputNode(WomIdentifier(parsedInputId), wdlTypeForInputParameter(input).get)
    }.toSet

    val workflowInputs: Map[String, GraphNodeOutputPort] =
      graphFromInputs.map {
        workflowInput =>
          workflowInput.localName -> workflowInput.singleOutputPort
      }.toMap

    val graphFromSteps: Checked[Set[GraphNode]] =
      steps.
        toList.
        foldLeft((Set.empty[GraphNode] ++ graphFromInputs).asRight[NonEmptyList[String]])(
          (nodes, step) => nodes.flatMap(step.callWithInputs(typeMap,  this, _, workflowInputs)))

    val graphFromOutputs: Checked[Set[GraphNode]] =
      outputs.toList.traverse[ErrorOr, GraphNode] {
        output =>

          val wdlType = cwlTypeToWdlType(output.`type`.flatMap(_.select[CwlType]).get)

          def lookupOutputSource(outputId: WorkflowOutputId): Checked[OutputPort] =
            for {
              set <- graphFromSteps
              call <- set.collectFirst { case callNode: CallNode if callNode.localName == outputId.stepId => callNode }.
                toRight(NonEmptyList.one(s"Call Node by name ${outputId.stepId} was not found in set $set"))
              output <- call.outputPorts.find(_.name == outputId.outputId).
                          toRight(NonEmptyList.one(s"looking for ${outputId.outputId} in call $call output ports ${call.outputPorts}"))
            } yield output

          lookupOutputSource(WorkflowOutputId(output.outputSource.flatMap(_.select[String]).get)).
            map(PortBasedGraphOutputNode(WomIdentifier(output.id), wdlType, _)).toValidated
      }.map(_.toSet).toEither

    for {
      outputs <- graphFromOutputs
      steps <- graphFromSteps
      ret <- Graph.validateAndConstruct(steps ++ graphFromInputs ++ outputs).toEither
    } yield ret
  }

  def womDefinition: Checked[WorkflowDefinition] = {
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
