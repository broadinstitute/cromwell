package wdl4s.cwl

import cats.data.NonEmptyList
import cats.syntax.either._
import shapeless._
import wdl4s.cwl.ScatterMethod._
import wdl4s.cwl.WorkflowStep.{Outputs, Run}
import wdl4s.wom.graph.GraphNodePort.{GraphNodeOutputPort, OutputPort}
import wdl4s.wom.graph.{CallNode, GraphNode, TaskCallNode}
import lenthall.Checked

import scala.language.postfixOps
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
                         hints: Option[Array[CwlAny]] = None,
                         label: Option[String] = None,
                         doc: Option[String] = None,
                         scatter: Option[String :+: Array[String] :+: CNil] = None,
                         scatterMethod: Option[ScatterMethod] = None) {


  def typedOutputs: WdlTypeMap = run.fold(RunOutputsToTypeMap)

  def fileName: Option[String] = run.select[String]

  val unqualifiedStepId = Try(WorkflowStepId(id)).map(_.stepId).getOrElse(id)

  /**
    * Generates all GraphNodes necessary to represent call nodes and input nodes
    * Recursive because dependencies are discovered as we iterate through inputs and corresponding
    * upstream nodes need to be generated on the fly.
    */
  def callWithInputs(typeMap: WdlTypeMap,
                     workflow: Workflow,
                     knownNodes: Set[GraphNode],
                     workflowInputs: Map[String, GraphNodeOutputPort]): Checked[Set[GraphNode]] = {

    // To avoid duplicating nodes, return immediately if we've already covered this node
    val haveWeSeenThisStep: Boolean = knownNodes.collect { case TaskCallNode(name, _, _, _) => name }.contains(unqualifiedStepId)

    if (haveWeSeenThisStep) Right(knownNodes)
    else {
      // Create a task definition for the underlying run.
      // For sub workflows, we'll need to handle the case where this could be a workflow definition
      //TODO: turn this select into a fold that supports other types of runnables
      val taskDefinition = run.select[CommandLineTool].map { _.taskDefinition } get

      /**
        * Method used to fold over the list of inputs declared by this step.
        * Note that because we work on saladed CWL, all ids are fully qualified at this point (e.g: file:///path/to/file/r.cwl#cgrep/pattern
        * The goal of this method is two fold (pardon the pun):
        *   1) link each input of the step to an output port (which at this point can be from a different step or from a workflow input)
        *   2) accumulate the nodes created along the way to achieve 1)
        */
      def foldInputs(mapAndNodes: Checked[(Map[String, OutputPort],  Set[GraphNode])],
                     workflowStepInput: WorkflowStepInput): Checked[(Map[String, OutputPort], Set[GraphNode])] =
        mapAndNodes flatMap {
        case (map, knownNodes) => //shadowing knownNodes on purpose to avoid accidentally referencing the outer one

          // The source from which we expect to satisfy this input (output from other step or workflow input)
          // TODO: this can be None in which case we should look if the "default" field is defined
          val inputSource: String = workflowStepInput.source.flatMap(_.select[String]).get

          /*
            * Try to find in the given set an output port named stepOutputId in a call node named stepId
            * This is useful when we've determined that the input points to an output of a different step and we want
            * to get the corresponding output port.
           */
          def findThisInputInSet(set: Set[GraphNode], stepId: String, stepOutputId: String): Checked[OutputPort] = {
            for {
            // We only care for outputPorts of call nodes
              call <- set.collectFirst { case callNode: CallNode if callNode.name == stepId => callNode }.
                        toRight(NonEmptyList.one(s"stepId $stepId not found in known Nodes $set"))
              output <- call.outputPorts.find(_.name == stepOutputId).
                          toRight(NonEmptyList.one(s"step output id $stepOutputId not found in ${call.outputPorts}"))
            } yield output
          }

          /*
            * Build a wom node for the given step and return the newly created nodes
            * This is useful when we've determined that the input belongs to an upstream step that we haven't covered yet
           */
          def buildUpstreamNodes(upstreamStepId: String): Checked[Set[GraphNode]] =
            // Find the step corresponding to this upstreamStepId in the set of all the steps of this workflow
            for {
              step <- workflow.steps.find { step => WorkflowStepId(step.id).stepId == upstreamStepId }.
                        toRight(NonEmptyList.one(s"no step of id $upstreamStepId found in ${workflow.steps.map(_.id)}"))
              call <- step.callWithInputs(typeMap, workflow, knownNodes, workflowInputs)
            } yield call

          // Parse the workflowStepInput to be able to access its components separately
          val parsedStepInputId = WorkflowStepInputOrOutputId(workflowStepInput.id)
          /*
            * Attempt to find the task definition input that corresponds to the workflow step input
            * Note that the Ids might be largely different (including the file name part)
            * We look for a task definition input that has the same name.
            * We assume here that workflow step inputs map 1 to 1 with definition (run) inputs and that they have the same name (last bit of the id)
            * e.g:
            * workflow step input could be: file:///Users/danb/wdl4s/r.cwl#cgrep/pattern
            * with a task description input: file:///Users/danb/wdl4s/r.cwl#cgrep/09f8bcac-a91a-49d5-afb6-2f1b1294e875/pattern
            *
            * The file name part could even be different if they have been saladed at a different time / place
           */
          val taskDefinitionInput = taskDefinition.inputs.map(_.name).find(_ == parsedStepInputId.ioId).getOrElse {
            throw new Exception(s"Can't find a corresponding input in the run section of ${parsedStepInputId.stepId} for ${parsedStepInputId.ioId}")
          }

          /*
            * Parse the inputSource (what this input is pointing to)
            * 2 cases:
            *   - points to a workflow input
            *   - points to an upstream step
           */
          FullyQualifiedName(inputSource) match {
            // The source points to a workflow input, which means it should be in the workflowInputs map
            case _: WorkflowInputId =>

              // Try to find it in the workflow inputs map, if we can't it's an error
              // TODO: (at least for now, it's possible to provide a default value but we don't support it yet)
              val outputPortMapping: Checked[(String, OutputPort)] = workflowInputs.collectFirst {
                case (inputId, port) if inputSource == inputId => taskDefinitionInput -> port
              }.toRight(NonEmptyList.one(s"Can't find workflow input for $inputSource"))

              outputPortMapping.map(mapping => (map + mapping) -> knownNodes)

            // The source points to an output from a different step
            case WorkflowStepInputOrOutputId(_, stepId, stepOutputId) =>
              val newNodesAndOutputPort: Checked[(Set[GraphNode], OutputPort)] =
              // First check if we've already built the WOM node for this step, and if so return the associated output port
                findThisInputInSet(knownNodes, stepId, stepOutputId).map(Set.empty[GraphNode] -> _)
                  .orElse {
                    // Otherwise build the upstream nodes
                    val newNodes = buildUpstreamNodes(stepId)
                    // And look again in those newly created nodes
                    newNodes.flatMap(nodes => findThisInputInSet(nodes, stepId, stepOutputId) map {
                      nodes -> _
                    })
                  }

              // Note that the key is the taskDefinitionInput (which in our previous example would be "pattern")
              // This is because WOM is going to use this map to link input ports to output ports
              // and will use the taskDefinition input port names to do so.
              // So we need to make sure that the keys in this map match the task definition input names (and not the workflow step input names !)
              newNodesAndOutputPort.map {
                case (newNodes, outputPort) => (map + (taskDefinitionInput -> outputPort), knownNodes ++ newNodes)
              }
          }
      }

      // We fold over inputs for this step and build an input -> output port lookup map as well as nodes created along the way
      val workflowOutputsMap: Checked[(Map[String, OutputPort], Set[GraphNode])] =
        in.foldLeft((Map.empty[String, OutputPort] -> knownNodes).asRight[NonEmptyList[String]]) (foldInputs)

      // Use what we've got to generate a call node and required input nodes
      // However here we expect WOM NOT to return any required input nodes because it would mean that some task definition inputs
      // have not been linked to either a workflow input or an upstream output, in which case they have no other way to be satisfied ( <- is that true ?)
      for {
        mapsAndNodes <- workflowOutputsMap
        (map, nodes) = mapsAndNodes
        call <-  CallNode.callWithInputs(unqualifiedStepId, taskDefinition, workflowInputs ++ map, Set.empty, prefixSeparator = "#").toEither
      } yield  call.nodes ++ nodes
    }
  }
}

/**
  * @see <a href="http://www.commonwl.org/v1.0/Workflow.html#WorkflowStepOutput">WorkflowstepOutput</a>
  */
case class WorkflowStepOutput(id: String)

object WorkflowStep {

  val emptyOutputs: Outputs = Coproduct[Outputs](Array.empty[String])

  type Run =
    String :+:
      CommandLineTool :+:
      ExpressionTool :+:
      Workflow :+:
      CNil

  type Outputs =
    Array[String] :+:
      Array[WorkflowStepOutput] :+:
      CNil
}
