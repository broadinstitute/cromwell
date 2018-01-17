package cwl

import java.nio.file.Paths

import cats.data.NonEmptyList
import cats.instances.list._
import cats.syntax.either._
import cats.syntax.traverse._
import common.Checked
import common.validation.Validation.OptionValidation
import common.validation.ErrorOr._
import common.validation.Validation._
import shapeless._
import shapeless.syntax.singleton._
import CwlVersion._
import cwl.command.ParentName
import wom.callable.WorkflowDefinition
import wom.executable.Executable
import wom.expression.{ValueAsAnExpression, WomExpression}
import wom.graph.GraphNodePort.{GraphNodeOutputPort, OutputPort}
import wom.graph._
import wom.types.WomType

case class Workflow private(
                             cwlVersion: Option[CwlVersion],
                             `class`: Witness.`"Workflow"`.T,
                             id: String,
                             inputs: Array[WorkflowInputParameter],
                             outputs: Array[WorkflowOutputParameter],
                             steps: Array[WorkflowStep],
                             requirements: Option[Array[Requirement]],
                             hints: Option[Array[Hint]]) {

  steps.foreach { _.parentWorkflow = this }

  /** Builds an `Executable` from a `Workflow` CWL with no parent `Workflow` */
  def womExecutable(validator: RequirementsValidator, inputFile: Option[String] = None): Checked[Executable] = {
    CwlExecutableValidation.buildWomExecutableCallable(womDefinition(validator, Vector.empty), inputFile)
  }

  // Circe can't create bidirectional links between workflow steps and runs (including `Workflow`s) so this
  // ugly var is here to link back to a possible parent workflow step. This is needed to navigate upward for finding
  // requirements in the containment hierarchy. There isn't always a containing workflow step so this is an `Option`.
  private[cwl] var parentWorkflowStep: Option[WorkflowStep] = None

  val allRequirements: List[Requirement] = requirements.toList.flatten ++ parentWorkflowStep.toList.flatMap { _.allRequirements }

  val allHints: List[Requirement] = {
    // Just ignore any hint that isn't a Requirement.
    val requirementHints = hints.toList.flatten.flatMap { _.select[Requirement] }
    requirementHints ++ parentWorkflowStep.toList.flatMap { _.allHints }
  }

  private [cwl] implicit val explicitWorkflowName = ParentName(id)

  val fileNames: List[String] = steps.toList.flatMap(_.run.select[String].toList)

  def outputsTypeMap: WomTypeMap = steps.foldLeft(Map.empty[String, WomType]) {
    // Not implemented as a `foldMap` because there is no semigroup instance for `WomType`s.  `foldMap` doesn't know that
    // we don't need a semigroup instance since the map keys should be unique and therefore map values would never need
    // to be combined under the same key.
    (acc, s) => acc ++ s.typedOutputs
  }

  def womGraph(workflowName: String, validator: RequirementsValidator, expressionLib: ExpressionLib): Checked[Graph] = {
    val workflowNameIdentifier = explicitWorkflowName.value.map(WomIdentifier.apply).getOrElse(WomIdentifier(workflowName))

    def womTypeForInputParameter(input: InputParameter): Option[WomType] = {
      input.`type`.map(_.fold(MyriadInputTypeToWomType))
    }

    val typeMap: WomTypeMap =
      outputsTypeMap ++
        // Note this is only looking at the workflow inputs and not recursing into steps, because our current thinking
        // is that in CWL graph inputs can only be defined at the workflow level.  It's possible that's not actually
        // correct, but that's the assumption being made here.
        inputs.toList.flatMap { i =>
          womTypeForInputParameter(i).map(i.id -> _).toList
        }.toMap

    val graphFromInputs: Set[ExternalGraphInputNode] = inputs.map {
      case inputParameter if inputParameter.default.isDefined =>
        val parsedInputId = FileAndId(inputParameter.id).id
        val womType = womTypeForInputParameter(inputParameter).get

        val womValue = inputParameter.default.get.fold(InputParameter.DefaultToWomValuePoly).apply(womType).toTry.get
        OptionalGraphInputNodeWithDefault(WomIdentifier(parsedInputId, inputParameter.id), womType, ValueAsAnExpression(womValue), parsedInputId)
      case input =>
        val parsedInputId = FileAndId(input.id).id

        RequiredGraphInputNode(WomIdentifier(parsedInputId, input.id), womTypeForInputParameter(input).get, parsedInputId)
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
          (nodes, step) => nodes.flatMap(step.callWithInputs(typeMap, this, _, workflowInputs, validator, expressionLib)))

    val graphFromOutputs: Checked[Set[GraphNode]] =
      outputs.toList.traverse[ErrorOr, GraphNode] {
        case WorkflowOutputParameter(id, _, _, _, _, _, _, Some(Inl(outputSource: String)), _, Some(tpe)) =>
          val womType:WomType = tpe.fold(MyriadOutputTypeToWomType)

          val parsedWorkflowOutput = FileAndId(id)
          val parsedOutputSource = FileStepAndId(outputSource)

          // Try to find an output port for this cwl output in the set of available nodes
          def lookupOutputSource(outputId: FileStepAndId): Checked[OutputPort] = {
            def isRightOutputPort(op: GraphNodePort.OutputPort) = FullyQualifiedName.maybeApply(op.name) match {
              case Some(fqn) => fqn.id == outputId.id
              case None => op.name == outputId.id
            }

            for {
              set <- graphFromSteps
              node <- set.collectFirst({
                case callNode: CallNode if callNode.localName == outputId.stepId => callNode
                case scatterNode: ScatterNode if scatterNode.innerGraph.calls.exists(_.localName == outputId.stepId) => scatterNode
              }).
                toRight(NonEmptyList.one(s"Call Node by name ${outputId.stepId} was not found in set $set"))
              output <- node.outputPorts.find(isRightOutputPort).toChecked(s"looking for ${outputId.id} in call $node output ports ${node.outputPorts}")
            } yield output
          }

          lookupOutputSource(parsedOutputSource).map({ port =>
            val localName = LocalName(parsedWorkflowOutput.id)
            val fullyQualifiedName = workflowNameIdentifier.fullyQualifiedName.combine(parsedWorkflowOutput.id)
            val outputIdentifier = WomIdentifier(localName, fullyQualifiedName)
            PortBasedGraphOutputNode(outputIdentifier, womType, port)
          }).toValidated
        case wop => throw new NotImplementedError(s"Workflow output parameters such as $wop are not supported.")
      }.map(_.toSet).toEither

    for {
      outputs <- graphFromOutputs
      steps <- graphFromSteps
      ret <- Graph.validateAndConstruct(steps ++ graphFromInputs ++ outputs).toEither
    } yield ret
  }

  def womDefinition(validator: RequirementsValidator, expressionLib: ExpressionLib): Checked[WorkflowDefinition] = {
    val name: String = Paths.get(id).getFileName.toString
    val meta: Map[String, String] = Map.empty
    val paramMeta: Map[String, String] = Map.empty
    val declarations: List[(String, WomExpression)] = List.empty

    womGraph(name, validator, expressionLib).map(graph =>
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

case class WorkflowInputParameter(id: String,
                                  label: Option[String] = None,
                                  secondaryFiles: Option[SecondaryFiles] = None,
                                  format: Option[Expression :+: String :+: Array[String] :+: CNil] = None,
                                  streamable: Option[Boolean] = None,
                                  doc: Option[String :+: Array[String] :+: CNil] = None,
                                  inputBinding: Option[InputCommandLineBinding] = None,
                                  default: Option[CwlAny] = None,
                                  `type`: Option[MyriadInputType] = None) extends InputParameter

object Workflow {

  def apply(cwlVersion: Option[CwlVersion] = Option(CwlVersion.Version1),
            id: String,
            inputs: Array[WorkflowInputParameter] = Array.empty,
            outputs: Array[WorkflowOutputParameter] = Array.empty,
            steps: Array[WorkflowStep] = Array.empty,
            requirements: Option[Array[Requirement]] = None,
            hints: Option[Array[Hint]] = None): Workflow  =
    Workflow(cwlVersion, "Workflow".narrow, id, inputs, outputs, steps, requirements, hints)
}
