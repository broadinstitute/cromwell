package cwl

import cats.Monoid
import cats.data.NonEmptyList
import cats.data.Validated._
import cats.instances.list._
import cats.syntax.either._
import cats.syntax.foldable._
import cats.syntax.monoid._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import cwl.ScatterLogic.ScatterVariablesPoly
import cwl.ScatterMethod._
import cwl.WorkflowStep.{WorkflowStepInputFold, _}
import cwl.WorkflowStepInput._
import cwl.command.ParentName
import shapeless.{:+:, CNil, _}
import wom.callable.Callable
import wom.callable.Callable._
import wom.graph.CallNode._
import wom.graph.GraphNodePort.{GraphNodeOutputPort, OutputPort}
import wom.graph._
import wom.graph.expression.ExpressionNode
import wom.types.WomType
import wom.values.WomValue

/**
  * An individual job to run.
  *
  * @see <a href="http://www.commonwl.org/v1.0/Workflow.html#WorkflowStep">CWL Spec | Workflow Step</a>
  * @param run Purposefully not defaulted as it's required in the specification and it is unreasonable to not have something to run.
  */
case class WorkflowStep(
                         id: String,
                         in: Array[WorkflowStepInput] = Array.empty,
                         out: WorkflowStepOutputType,
                         run: Run,
                         requirements: Option[Array[Requirement]] = None,
                         hints: Option[Array[Hint]] = None,
                         label: Option[String] = None,
                         doc: Option[String] = None,
                         scatter: ScatterVariables = None,
                         scatterMethod: Option[ScatterMethod] = None) {

  run.select[Workflow].foreach(_.parentWorkflowStep = Option(this))
  run.select[CommandLineTool].foreach(_.parentWorkflowStep = Option(this))
  run.select[ExpressionTool].foreach(_.parentWorkflowStep = Option(this))

  // We're scattering if scatter is defined, and if it's a list of variables the list needs to be non empty
  private val isScattered: Boolean = scatter exists { _.select[Array[String]].forall(_.nonEmpty) }

  // Circe can't create bidirectional links between workflows and workflow steps so this ugly var is here to link
  // back to the parent workflow. This is needed to navigate upward for finding requirements in the containment
  // hierarchy. There is always a workflow containing a workflow step so this is not an `Option`.
  private[cwl] var parentWorkflow: Workflow = _

  lazy val allRequirements = RequirementsAndHints(requirements.toList.flatten ++ parentWorkflow.allRequirements.list)

  lazy val womFqn: wom.graph.FullyQualifiedName = {
    implicit val parentName = parentWorkflow.explicitWorkflowName
    val localFqn = FullyQualifiedName.maybeApply(id).map(_.id).getOrElse(id)
    parentWorkflow.womFqn.map(_.combine(localFqn)).getOrElse(wom.graph.FullyQualifiedName(localFqn))
  }

  lazy val allHints: List[Requirement] = {
    // Just ignore any hint that isn't a Requirement.
    val requirementHints = hints.toList.flatten.flatMap { _.select[Requirement] }
    requirementHints ++ parentWorkflow.allHints
  }

  // If the step is being scattered over, apply the necessary transformation to get the final output array type.
  lazy val scatterTypeFunction: WomType => WomType = scatter.map(_.fold(ScatterVariablesPoly)) match {
    case Some(Nil) => identity[WomType]
    case Some(nonEmpty) => ScatterLogic.scatterGatherPortTypeFunction(scatterMethod, NonEmptyList.fromListUnsafe(nonEmpty))
    case _ =>  identity[WomType]
  }

  def typedOutputs: WomTypeMap = {
    implicit val parentName = ParentName(id)
    // Find the type of the outputs of the run section
    val runOutputTypes = run.fold(RunOutputsToTypeMap).apply(allRequirements.schemaDefRequirement)
      .map({
        case (runOutputId, womType) => FullyQualifiedName(runOutputId).id -> womType
      })
    // Use them to find get the final type of the workflow outputs, and only the workflow outputs
    out.map({ stepOutput =>
      val stepOutputValue = stepOutput.select[WorkflowStepOutput].map(_.id).getOrElse(stepOutput.select[String].get)
      val stepOutputId = FullyQualifiedName(stepOutputValue)
      stepOutputValue -> scatterTypeFunction(runOutputTypes(stepOutputId.id))
    }).toMap
  }

  def fileName: Option[String] = run.select[String]

  /**
    * Generates all GraphNodes necessary to represent call nodes and input nodes
    * Recursive because dependencies are discovered as we iterate through inputs and corresponding
    * upstream nodes need to be generated on the fly.
    *
    * Example:
    *
    * CWL:
    * inputs:
    *    -id: workflow_input_A 
    *    type: string[]
    *   -id: workflow_input_B 
    *    type: string[]
    *    -id: workflow_input_C 
    *    type: string  
    * steps: 
    *   -id: echo
    *    run: echo.cwl
    *    scatter: input0
    *    in:
    *      -id:input0
    *       source:
    *         - "#workflow_input_A"
    *         - "#workflow_input_B" 
    *       valueFrom:"bonjour"
    *
    *      -id: input1
    *        source: "#workflow_input_C"
    *       valueFrom:"$(inputs.input0)"
    *
    * WOM Diagram:
    *
    *  +------------+   +------------+   +------------+
    *  | WF Input A |   | WF Input B |   | WF Input C |
    *  +--------+---+ +-+------------+   +------+-----+
    *           |     |                         |
    *           |     |                         |
    *     +-----v-----v+                 +------v-----+
    *     | StepInput0 |                 | StepInput1 |
    *     | MergeLogic |                 | MergeLogic |
    *     +-----+------+                 +----+-------+
    *           |                             |
    *  +----------------------------------------------+
    *  |        |                             |       |
    *  | +------v----------+              +---v--+    |
    *  | | ScatterVariable +------------+ | OGIN |    |
    *  | +-----+-----------+            | +--+-+-+    |
    *  |       |                        |    | |      |
    *  | +-----v----+                   |    | |      |
    *  | |StepInput0<------------------------+ |      |
    *  | |Expression|                   |      |      |
    *  | +----+-----+                +--v------v+     |
    *  |      |                      |StepInput1|     |
    *  |      |                      |Expression|     |
    *  |      |                      +------+---+     |
    *  |      |    +-----------+            |         |
    *  |      +----> Call Node <------------+         |
    *  |           +-----------+                      |
    *  |                                 Scatter Node |
    *  +----------------------------------------------+
    *
    * MergeNode: If the step input has one or more sources, a merge node will be created and responsible for merging
    *   those input sources together. It will NOT evaluate the valueFrom field of the input.
    *
    * ScatterVariableNode: If the step input is being scattered over, a scatter variable node will be created and will
    *   act as a proxy inside the scatter graph for the shards of the scatter. It depends on an upstream merge node outputing an array
    *   and will provide at runtime shard values for other nodes of the scatter graph.
    *
    * OGIN: If the step has at least one input being scattered over, there will be a scatter node created.
    *   For inputs that are NOT being scattered over but still have one or more input sources (and hence a merge node), an OGIN
    *   will be created to act as a proxy to the merge node outside the scatter graph.
    *
    * ExpressionNode: If an input has a valueFrom field, an expression node will be created to evaluate the expression.
    *   An important fact to note is that the expression needs access to all other input values
    *   AFTER their source, default value and shard number has been determined but
    *   BEFORE their (potential) valueFrom is evaluated (see http://www.commonwl.org/v1.0/Workflow.html#WorkflowStepInput)
    *   This is why on the above diagram, StepInput0Expression depends on the OGIN, and StepInput1Expression depends on the scatter variable.
    */
  def callWithInputs(typeMap: WomTypeMap,
                     workflow: Workflow,
                     knownNodes: Set[GraphNode],
                     workflowInputs: Map[String, GraphNodeOutputPort],
                     validator: RequirementsValidator,
                     expressionLib: ExpressionLib): Checked[Set[GraphNode]] = {

    implicit val parentName = workflow.explicitWorkflowName

    val scatterLookupSet =
      scatter.toList.
        flatMap(_.fold(StringOrStringArrayToStringList)).
        map(id => FullyQualifiedName(id).id)

    def isStepScattered(workflowStepInputId: String) = scatterLookupSet.contains(workflowStepInputId)

    val unqualifiedStepId: WomIdentifier = {
      FullyQualifiedName.maybeApply(id).map({ fqn =>
        WomIdentifier(LocalName(fqn.id), womFqn)
      }).getOrElse(WomIdentifier(id))
    }

    def typedRunInputs: Map[String, Option[MyriadInputType]] = run.fold(RunToInputTypeMap).apply(parentName)

    def allIdentifiersRecursively(nodes: Set[GraphNode]): Set[WomIdentifier] = nodes.flatMap({
      case w: WorkflowCallNode=> Set(w.identifier)
      case c: CommandCallNode => Set(c.identifier)
      case e: ExpressionCallNode => Set(e.identifier)
      // When a node a call node is being scattered over, it is wrapped inside a scatter node. We still don't want to
      // duplicate it though so look inside scatter nodes to see if it's there.
      case scatter: ScatterNode => allIdentifiersRecursively(scatter.innerGraph.nodes)
      case _ => Set.empty[WomIdentifier]
    })

    // To avoid duplicating nodes, return immediately if we've already covered this node
    val haveWeSeenThisStep: Boolean = allIdentifiersRecursively(knownNodes).contains(unqualifiedStepId)

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
          /*
            * Try to find in the given set an output port named stepOutputId in a call node named stepId
            * This is useful when we've determined that the input points to an output of a different step and we want
            * to get the corresponding output port.
           */
          def findThisInputInSet(set: Set[GraphNode], stepId: String, stepOutputId: String): Checked[OutputPort] = {
            for {
              // We only care for outputPorts of call nodes or scatter nodes
              call <- set.collectFirst {
                case callNode: CallNode if callNode.localName == stepId => callNode
                case scatterNode: ScatterNode if scatterNode.innerGraph.calls.exists(_.localName == stepId) => scatterNode
              }.
                toRight(NonEmptyList.one(s"stepId $stepId not found in known Nodes $set"))
              output <- call.outputPorts.find(_.internalName == stepOutputId).
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
            findThisInputInSet(accumulatedNodes, stepId, stepOutputId).map(outputPort => (Map(s"$stepId/$stepOutputId" -> outputPort), accumulatedNodes))
              .orElse {
                // Otherwise build the upstream nodes and look again in those newly created nodes
                for {
                  newNodes <- buildUpstreamNodes(stepId, accumulatedNodes)
                  sourceMappings <- findThisInputInSet(newNodes, stepId, stepOutputId).map(outputPort => Map(s"$stepId/$stepOutputId" -> outputPort))
                } yield (sourceMappings, newNodes ++ accumulatedNodes)
              }
          }

          lazy val workflowStepInputId = FullyQualifiedName(workflowStepInput.id).id

          def updateFold(sourceMappings: Map[String, OutputPort], newNodes: Set[GraphNode]): Checked[WorkflowStepInputFold] = {
            val typeExpectedByRunInput: Option[cwl.MyriadInputType] = typedRunInputs.get(workflowStepInputId).flatten

            val isThisStepScattered = isStepScattered(workflowStepInputId)

            workflowStepInput.toMergeNode(sourceMappings, expressionLib, typeExpectedByRunInput, isThisStepScattered, allRequirements.schemaDefRequirement) match {
              // If the input needs a merge node, build it and add it to the input fold
              case Some(mergeNode) =>
                mergeNode.toEither.map({ node =>
                  fold |+| WorkflowStepInputFold(
                    mergeNodes = Map(workflowStepInput -> node),
                    generatedNodes = newNodes
                  )
                })
              case None => (fold |+| WorkflowStepInputFold(generatedNodes = newNodes)).validNelCheck
            }
          }

          /*
           * We intend to validate that all of these sources point to a WOM Outputport that we know about.
           *
           * If we don't know about them, we find upstream nodes and build them (see "buildUpstreamNodes").
           */
          val baseCase = (Map.empty[String, OutputPort], fold.generatedNodes).asRight[NonEmptyList[String]]
          val inputMappingsAndGraphNodes: Checked[(Map[String, OutputPort], Set[GraphNode])] =
            workflowStepInput.sources.foldLeft(baseCase) {
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
                  case FileStepAndId(_, _, stepId, stepOutputId) => fromStepOutput(stepId, stepOutputId, graphNodes).map({ case (newMap, newNodes) => (sourceMappings ++ newMap, newNodes) })
                }
              case (other, _) => other
            }

          inputMappingsAndGraphNodes.flatMap((updateFold _).tupled)
      }

      /*
       * Folds over input definitions and build an InputDefinitionFold
       */
      def foldInputDefinition(pointerNode: Map[String, GraphNodeWithSingleOutputPort])
                             (inputDefinition: InputDefinition): ErrorOr[InputDefinitionFold] = {
        inputDefinition match {
          case _ if pointerNode.contains(inputDefinition.name) =>
            val expressionNode = pointerNode(inputDefinition.name)
            InputDefinitionFold(
              mappings = List(inputDefinition -> expressionNode.inputDefinitionPointer),
              callInputPorts = Set(callNodeBuilder.makeInputPort(inputDefinition, expressionNode.singleOutputPort))
            ).validNel

          // No expression node mapping, use the default
          case withDefault @ OverridableInputDefinitionWithDefault(_, _, expression, _, _) =>
            InputDefinitionFold(
              mappings = List(withDefault -> Coproduct[InputDefinitionPointer](expression))
            ).validNel

          // Required input without default value and without mapping, this is a validation error
          case RequiredInputDefinition(requiredName, _, _, _) =>
            s"Input ${requiredName.value} is required and is not bound to any value".invalidNel

          // Optional input without mapping, defaults to empty value
          case optional: OptionalInputDefinition =>
            InputDefinitionFold(
              mappings = List(optional -> Coproduct[InputDefinitionPointer](optional.womType.none: WomValue))
            ).validNel
        }
      }

      /*
        If the step is being scattered over, then merge nodes can't directly be referenced because they will be outside the scatter graph.
        For inputs that are being scattered over, a scatter variable has already been created, but for the others we need
        an OGIN to link the merge node to the inner scatter graph.
       */
      def buildOGINs(mergeNodes: Map[WorkflowStepInput, ExpressionNode],
                     scatterVariables: Map[WorkflowStepInput, ScatterVariableNode]): Map[WorkflowStepInput, OuterGraphInputNode] = if (isScattered) {
        mergeNodes
          .collect({
            case (input, mergeNode) if !scatterVariables.contains(input) =>
              val ogin = OuterGraphInputNode(
                WomIdentifier(input.parsedId).combine("OGIN"),
                mergeNode.singleOutputPort,
                preserveScatterIndex = false
              )
              input -> ogin
          })
      } else Map.empty

      /*
        * For inputs that have a valueFrom field, create an ExpressionNode responsible for evaluating the expression.
        * Note that this expression might need access to the other input values, so make each expression node depend on all other
        * inputs.
       */
      def buildStepInputValueFromNodes(sharedInputNodes: Map[WorkflowStepInput, GraphNodeWithSingleOutputPort]): Checked[Map[String, ExpressionNode]] = {
        // Add new information to the typeMap from the shard input nodes.
        lazy val updatedTypeMap = sharedInputNodes.map({
          // If the input node is a scatter variable, make sure the type is the item type, not the array type, as the expression node
          // will operate on shards not on the whole scattered array.
          case (stepInput, scatter: ScatterVariableNode) => stepInput.parsedId -> scatter.womType
          case (stepInput, node) => stepInput.parsedId -> node.singleOutputPort.womType
        }) ++ typeMap

        // Go over each step input and create an expression node for those which have a valueFrom
        in.toList.collect({
          case stepInput @ WorkflowStepInput(_, _, _, _, Some(valueFrom)) =>
            // Transform the shared inputs map into a usable map to create the expression.
            lazy val sharedInputMap: Map[String, OutputPort] = sharedInputNodes.map({
              case (siblingStepInput, graphNode) => siblingStepInput.parsedId -> graphNode.singleOutputPort
            })
            val typeExpectedByRunInput: Option[cwl.MyriadInputType] = typedRunInputs.get(stepInput.parsedId).flatten
            val isThisStepScattered = isStepScattered(stepInput.parsedId)

            stepInput.toExpressionNode(valueFrom, typeExpectedByRunInput, isThisStepScattered, sharedInputMap, updatedTypeMap, expressionLib, allRequirements.schemaDefRequirement).map(stepInput.parsedId -> _)
        })
          .sequence[ErrorOr, (String, ExpressionNode)]
          .toEither
          .map(_.toMap)
      }

      //inputs base case consist of the nodes we already know about
      val baseCase = WorkflowStepInputFold(generatedNodes = knownNodes).asRight[NonEmptyList[String]]

      // WorkflowStepInputFold contains the mappings from step input to ExpressionNode as well as all created nodes
      val stepInputFoldCheck: Checked[WorkflowStepInputFold] = in.foldLeft(baseCase)(foldStepInput)

      /*
        * This (big) flatMap builds nodes from top to bottom in the diagram above.
        * If necessary, the scatter node is built last as it wraps some of the other nodes.
       */
      for {
        /* ************************************ */
        /* ************ Merge Nodes *********** */
        /* ************************************ */
        // Build merge nodes and recursively generates other call nodes that we haven't seen so far
        stepInputFold <- stepInputFoldCheck
        // Extract the merge nodes from the fold
        mergeNodes = stepInputFold.mergeNodes

        /* ************************************ */
        /* ****** Scatter Variable Nodes ****** */
        /* ************************************ */
        scatterVariableNodes <- ScatterLogic.buildScatterVariableNodes(scatter, mergeNodes, unqualifiedStepId.localName.value)

        /* ************************************ */
        /* *************** OGINS ************** */
        /* ************************************ */
        ogins = buildOGINs(mergeNodes, scatterVariableNodes)

        /* ************************************ */
        /* ********* Expression Nodes ********* */
        /* ************************************ */
        // Aggregate the generated nodes so far. This map will be used to generate expression nodes, so the order of aggregation matters:
        // scatter variables and ogins take precedence over merge nodes (see diagram)
        aggregatedMapForValueFromNodes = mergeNodes ++ scatterVariableNodes ++ ogins
        // Build expression nodes for inputs that have a valueFrom field
        stepInputValueFromNodes <- buildStepInputValueFromNodes(aggregatedMapForValueFromNodes)

        /* ************************************ */
        /* ************* Call Node ************ */
        /* ************************************ */
        // Get the callable object for this step
        checkedCallable <- callable
        // Aggregate again by adding generated expression nodes. Again order matters here, expression nodes override other nodes.
        aggregatedMapForInputDefinitions = aggregatedMapForValueFromNodes.asIdentifierMap ++ stepInputValueFromNodes
        // Assign each of the callable's input definition to an output port from the pointer map
        inputDefinitionFold <- checkedCallable.inputs.foldMap(foldInputDefinition(aggregatedMapForInputDefinitions)).toEither
        // Build the call node
        callAndNodes = callNodeBuilder.build(unqualifiedStepId, checkedCallable, inputDefinitionFold, Set.empty, None)
        // Depending on whether the step is being scattered, invoke the scatter node builder or not

        /* ************************************ */
        /* ************ Scatter Node ********** */
        /* ************************************ */
        scatterNodeOrExposedNodes <- if (isScattered) {
          ScatterLogic.buildScatterNode(
            callAndNodes,
            NonEmptyList.fromListUnsafe(scatterVariableNodes.values.toList),
            ogins.values.toSet,
            stepInputValueFromNodes.values.toSet,
            scatterMethod).map(Set(_))
        } else {
          // If there's no scatter node then we need to return the expression nodes and the call node explicitly
          // as they won't be contained in the scatter inner graph
          (stepInputValueFromNodes.values.toSet + callAndNodes.node).validNelCheck
        }

        /*
          * Return all the nodes that need to be made available to the workflow graph:
          * knownNodes: this method is used to fold over steps so we don't want to forget to accumulate known nodes
          * mergeNodes: they're always outside of the scatter so always return them
          * generatedNodes: nodes generated recursively to build this node
          * scatterNodeOrExposedNodes: see explanation above
         */
        allNodes = knownNodes ++ mergeNodes.values.toSet ++ stepInputFold.generatedNodes ++ scatterNodeOrExposedNodes
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
        mergeNodes = x.mergeNodes ++ y.mergeNodes,
        generatedNodes = x.generatedNodes ++ y.generatedNodes
      )
    }
  }

  private [cwl] case class WorkflowStepInputFold(mergeNodes: Map[WorkflowStepInput, ExpressionNode] = Map.empty,
                                                 generatedNodes: Set[GraphNode] = Set.empty)

  /**
    * Maps input variable (to be scattered over) to their scatter variable node
    */
  type ScatterMappings = Map[ExpressionNode, ScatterVariableNode]

  val emptyOutputs: WorkflowStepOutputType = Array.empty

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

  type WorkflowStepOutputInnerType = String :+: WorkflowStepOutput :+: CNil
  type WorkflowStepOutputType = Array[WorkflowStepOutputInnerType]
}
