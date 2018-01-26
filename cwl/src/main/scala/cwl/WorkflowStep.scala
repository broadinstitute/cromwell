package cwl

import cats.Monoid
import cats.data.NonEmptyList
import cats.data.Validated._
import cats.instances.list._
import cats.syntax.either._
import cats.syntax.foldable._
import cats.syntax.monoid._
import cats.syntax.validated._
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import cwl.ScatterMethod._
import cwl.WorkflowStep.{WorkflowStepInputFold, _}
import shapeless._
import wom.callable.Callable
import wom.callable.Callable._
import wom.graph.CallNode._
import wom.graph.GraphNodePort.{GraphNodeOutputPort, OutputPort}
import wom.graph._
import wom.graph.expression.ExpressionNode
import wom.values.WomValue

import scala.util.Try

/**
  * An individual job to run.
  *
  * @see <a href="http://www.commonwl.org/v1.0/Workflow.html#WorkflowStep">CWL Spec | Workflow Step</a>
  * @param run Purposefully not defaulted as it's required in the specification and it is unreasonable to not have something to run.
  */
case class WorkflowStep(
                         id: String,
                         in: Array[WorkflowStepInput] = Array.empty,
                         out: Outputs,
                         run: Run,
                         requirements: Option[Array[Requirement]] = None,
                         hints: Option[Array[Hint]] = None,
                         label: Option[String] = None,
                         doc: Option[String] = None,
                         scatter: ScatterVariables = None,
                         scatterMethod: Option[ScatterMethod] = None) {

  run.select[Workflow].foreach(_.parentWorkflowStep = Option(this))
  run.select[CommandLineTool].foreach(_.parentWorkflowStep = Option(this))

  // We're scattering if scatter is defined, and if it's a list of variables the list needs to be non empty
  private val isScattered: Boolean = scatter exists { _.select[Array[String]].forall(_.nonEmpty) }

  // Circe can't create bidirectional links between workflows and workflow steps so this ugly var is here to link
  // back to the parent workflow. This is needed to navigate upward for finding requirements in the containment
  // hierarchy. There is always a workflow containing a workflow step so this is not an `Option`.
  private[cwl] var parentWorkflow: Workflow = _

  lazy val allRequirements: List[Requirement] = requirements.toList.flatten ++ parentWorkflow.allRequirements

  lazy val allHints: List[Requirement] = {
    // Just ignore any hint that isn't a Requirement.
    val requirementHints = hints.toList.flatten.flatMap { _.select[Requirement] }
    requirementHints ++ parentWorkflow.allHints
  }

  def typedOutputs: WomTypeMap = run.fold(RunOutputsToTypeMap)

  def fileName: Option[String] = run.select[String]

  /**
    * Generates all GraphNodes necessary to represent call nodes and input nodes
    * Recursive because dependencies are discovered as we iterate through inputs and corresponding
    * upstream nodes need to be generated on the fly.
    */
  def callWithInputs(typeMap: WomTypeMap,
                     workflow: Workflow,
                     knownNodes: Set[GraphNode],
                     workflowInputs: Map[String, GraphNodeOutputPort],
                     validator: RequirementsValidator,
                     expressionLib: ExpressionLib): Checked[Set[GraphNode]] = {

    implicit val parentName = workflow.explicitWorkflowName
    
    val unqualifiedStepId = WomIdentifier(Try(FullyQualifiedName(id)).map(_.id).getOrElse(id))

    // To avoid duplicating nodes, return immediately if we've already covered this node
    val haveWeSeenThisStep: Boolean = knownNodes.collect {
      case CommandCallNode(identifier, _, _, _) => identifier
      case ExpressionCallNode(identifier, _, _, _) => identifier
      case WorkflowCallNode(identifier, _, _, _) => identifier
    }.contains(unqualifiedStepId)

    if (haveWeSeenThisStep) Right(knownNodes)
    else {
      val callable: Checked[Callable] = run match {
        case Run.CommandLineTool(clt) => clt.buildTaskDefinition(validator, expressionLib)
        case Run.Workflow(wf) => wf.womDefinition(validator, expressionLib)
        case Run.ExpressionTool(et) => et.buildTaskDefinition(validator, expressionLib)
      }

      val callNodeBuilder = new CallNode.CallNodeBuilder()

      /*
       * Method used to fold over the list of inputs declared by this step.
       * Note that because we work on saladed CWL, all ids are fully qualified at this point (e.g: file:///path/to/file/three_step.cwl#cgrep/pattern
       * The goal of this method is two fold (pardon the pun):
       *   1) link each input of the step to an output port (which at this point can be from a different step or from a workflow input)
       *   2) accumulate the nodes created along the way to achieve 1)
       */
      def foldStepInput(currentFold: Checked[WorkflowStepInputFold], workflowStepInput: WorkflowStepInput): Checked[WorkflowStepInputFold] = currentFold flatMap {
        fold =>
          // The sources from which we expect to satisfy this input (output from other step or workflow input)
          val inputSources: List[String] = workflowStepInput.source.toList.flatMap(_.fold(WorkflowStepInputSourceToStrings))

          /*
            * Try to find in the given set an output port named stepOutputId in a call node named stepId
            * This is useful when we've determined that the input points to an output of a different step and we want
            * to get the corresponding output port.
           */
          def findThisInputInSet(set: Set[GraphNode], stepId: String, stepOutputId: String): Checked[OutputPort] = {
            for {
              // We only care for outputPorts of call nodes
              call <- set.collectFirst { case callNode: CallNode if callNode.localName == stepId => callNode }.
                toRight(NonEmptyList.one(s"stepId $stepId not found in known Nodes $set"))
              output <- call.outputPorts.find(_.name == stepOutputId).
                toRight(NonEmptyList.one(s"step output id $stepOutputId not found in ${call.outputPorts}"))
            } yield output
          }

          /*
           * Build a wom node for the given step and return the newly created nodes
           * This is useful when we've determined that the input belongs to an upstream step that we haven't covered yet
           */
          def buildUpstreamNodes(upstreamStepId: String, accumulatedNodes: Set[GraphNode]): Checked[Set[GraphNode]] =
          // Find the step corresponding to this upstreamStepId in the set of all the steps of this workflow
            for {
              step <- workflow.steps.find { step => FullyQualifiedName(step.id).id == upstreamStepId }.
                toRight(NonEmptyList.one(s"no step of id $upstreamStepId found in ${workflow.steps.map(_.id).toList}"))
              call <- step.callWithInputs(typeMap, workflow, accumulatedNodes, workflowInputs, validator, expressionLib)
            } yield call

          def fromWorkflowInput(inputName: String): Checked[Map[String, OutputPort]] = {
            // Try to find it in the workflow inputs map, if we can't it's an error
            workflowInputs.collectFirst {
              case (inputId, port) if inputName == inputId => Map(inputId -> port).asRight[NonEmptyList[String]]
            } getOrElse s"Can't find workflow input for $inputName".invalidNelCheck[Map[String, OutputPort]]
          }

          def fromStepOutput(stepId: String, stepOutputId: String, accumulatedNodes: Set[GraphNode]): Checked[(Map[String, OutputPort], Set[GraphNode])] = {
            // First check if we've already built the WOM node for this step, and if so return the associated output port
            findThisInputInSet(accumulatedNodes, stepId, stepOutputId).map(outputPort => (Map(stepOutputId -> outputPort), Set.empty[GraphNode]))
              .orElse {
                // Otherwise build the upstream nodes and look again in those newly created nodes
                for {
                  newNodes <- buildUpstreamNodes(stepId, accumulatedNodes)
                  sourceMappings <- findThisInputInSet(newNodes, stepId, stepOutputId).map(outputPort => Map(stepOutputId -> outputPort))
                } yield (sourceMappings, newNodes)
              }
          }

          def updateFold(sourceMappings: Map[String, OutputPort], newNodes: Set[GraphNode]): Checked[WorkflowStepInputFold] =
            workflowStepInput.toExpressionNode(sourceMappings, typeMap, expressionLib).map({ expressionNode =>
              fold |+| WorkflowStepInputFold(
                stepInputMapping = Map(FullyQualifiedName(workflowStepInput.id).id -> expressionNode),
                generatedNodes = newNodes + expressionNode
              )
            }).toEither

          val baseCase = (Map.empty[String, OutputPort], fold.generatedNodes).asRight[NonEmptyList[String]]
          val inputMappingsAndGraphNodes: Checked[(Map[String, OutputPort], Set[GraphNode])] =
            inputSources.foldLeft(baseCase){
              case (Right((sourceMappings, graphNodes)), inputSource) =>
                /*
               * Parse the inputSource (what this input is pointing to)
               * 2 cases:
               *   - points to a workflow input
               *   - points to an upstream step
               */
                FullyQualifiedName(inputSource) match {
                  // The source points to a workflow input, which means it should be in the workflowInputs map
                  case FileAndId(_, _, inputId) => fromWorkflowInput(inputId).map(newMap => (sourceMappings ++ newMap, graphNodes))
                  // The source points to an output from a different step
                  case FileStepAndId(_, _, stepId, stepOutputId) => fromStepOutput(stepId, stepOutputId, graphNodes)
                }
              case (other, _) => other
            }

          inputMappingsAndGraphNodes.flatMap((updateFold _).tupled)
      }

      /*
       * Folds over input definitions and build an InputDefinitionFold
       */
      def foldInputDefinition(expressionNodes: Map[String, ExpressionNode], scatterVariables: List[ScatterVariableNode])
                             (inputDefinition: InputDefinition): ErrorOr[InputDefinitionFold] = {
        inputDefinition match {
          // First 2 cases: We got an expression node, meaning there was a workflow step input for this input definition
          // Depending on whether or not this input is being scattered over take appropriate action
          case _ if expressionNodes.contains(inputDefinition.name) && isScattered =>
            val expressionNode = expressionNodes(inputDefinition.name)
            ScatterLogic.buildScatteredInputFold(inputDefinition, scatterVariables, expressionNode, callNodeBuilder)

          case _ if expressionNodes.contains(inputDefinition.name) =>
            val expressionNode = expressionNodes(inputDefinition.name)
            InputDefinitionFold(
              mappings = List(inputDefinition -> expressionNode.inputDefinitionPointer),
              callInputPorts = Set(callNodeBuilder.makeInputPort(inputDefinition, expressionNode.singleExpressionOutputPort)),
              newExpressionNodes = Set(expressionNode)
            ).validNel

          // No expression node mapping, use the default
          case withDefault @ InputDefinitionWithDefault(_, _, expression, _) =>
            InputDefinitionFold(
              mappings = List(withDefault -> Coproduct[InputDefinitionPointer](expression))
            ).validNel

          // Required input without default value and without mapping, this is a validation error
          case RequiredInputDefinition(requiredName, _, _) =>
            s"Input $requiredName is required and is not bound to any value".invalidNel

          // Optional input without mapping, defaults to empty value
          case optional: OptionalInputDefinition =>
            InputDefinitionFold(
              mappings = List(optional -> Coproduct[InputDefinitionPointer](optional.womType.none: WomValue))
            ).validNel
        }
      }

      //inputs base case consist of the nodes we already know about
      val baseCase = WorkflowStepInputFold(generatedNodes = knownNodes).asRight[NonEmptyList[String]]

      // WorkflowStepInputFold contains the mappings from step input to ExpressionNode as well as all created nodes
      val stepInputFoldCheck: Checked[WorkflowStepInputFold] = in.foldLeft(baseCase)(foldStepInput)

      /*
        1) Fold over the workflow step inputs:
          - Create an expression node for each input
          - recursively generates unseen call nodes as we discover them going through step input sources
          - accumulate all that in the WorkflowStepInputFold
        2) Fold over the callable input definition using the expression node map from 1):
          - determine the correct mapping for the input definition based on the expression node map
          and the type of input definition
          - accumulate those mappings, along with potentially newly created graph input nodes as well as call input ports
          in an InputDefinitionFold
        3) Use the InputDefinitionFold to build a new call node
       */
      for {
        stepInputFold <- stepInputFoldCheck
        checkedCallable <- callable
        scatterVariables <- ScatterLogic.buildScatterVariables(scatter, stepInputFold, unqualifiedStepId.localName.value)
        inputDefinitionFold <- checkedCallable.inputs.foldMap(foldInputDefinition(stepInputFold.stepInputMapping, scatterVariables)).toEither
        callAndNodes = callNodeBuilder.build(unqualifiedStepId, checkedCallable, inputDefinitionFold)
        allNodes <- if (scatterVariables.nonEmpty) {
          // NB: Pattern matching the scatterVariables against "head :: tail" to build a NeL doesn't work because the compiler takes :: for the shapeless operator
          ScatterLogic.prepareNodesForScatteredCall(knownNodes, callAndNodes, stepInputFold, NonEmptyList.fromListUnsafe(scatterVariables), scatterMethod)
        } else {
          (knownNodes ++ callAndNodes.nodes ++ stepInputFold.generatedNodes).validNelCheck
        }
      } yield allNodes
    }
  }
}

/**
  * @see <a href="http://www.commonwl.org/v1.0/Workflow.html#WorkflowStepOutput">WorkflowstepOutput</a>
  */
case class WorkflowStepOutput(id: String)

object WorkflowStep {

  // A monoid can't be derived automatically for this class because it contains a Map[String, ExpressionNode],
  // and there's no monoid defined over ExpressionNode
  implicit val workflowStepInputFoldMonoid: Monoid[WorkflowStepInputFold] = new Monoid[WorkflowStepInputFold] {
    override def empty: WorkflowStepInputFold = WorkflowStepInputFold()
    override def combine(x: WorkflowStepInputFold, y: WorkflowStepInputFold): WorkflowStepInputFold = {
      WorkflowStepInputFold(
        stepInputMapping = x.stepInputMapping ++ y.stepInputMapping,
        generatedNodes = x.generatedNodes ++ y.generatedNodes
      )
    }
  }

  private [cwl] case class WorkflowStepInputFold(stepInputMapping: Map[String, ExpressionNode] = Map.empty,
                                                 generatedNodes: Set[GraphNode] = Set.empty)

  /**
    * Maps input variable (to be scattered over) to their scatter variable node
    */
  type ScatterMappings = Map[ExpressionNode, ScatterVariableNode]

  val emptyOutputs: Outputs = Coproduct[Outputs](Array.empty[String])

  type Run =
    String :+:
      CommandLineTool :+:
      ExpressionTool :+:
      Workflow :+:
      CNil

  object Run {
    object String { def unapply(run: Run): Option[String] = run.select[String] }
    object Workflow { def unapply(run: Run): Option[Workflow] = run.select[Workflow] }
    object CommandLineTool { def unapply(run: Run): Option[CommandLineTool] = run.select[CommandLineTool] }
    object ExpressionTool { def unapply(run: Run): Option[ExpressionTool] = run.select[ExpressionTool] }
  }

  type Outputs =
    Array[String] :+:
      Array[WorkflowStepOutput] :+:
      CNil
}
