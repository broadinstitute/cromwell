package wdl4s.wom.graph

import cats.syntax.traverse._
import cats.instances.list._
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.wom.callable.Callable.OutputDefinition
import wdl4s.wom.callable.{Callable, TaskDefinition, WorkflowDefinition}
import wdl4s.wom.graph.GraphNode.{GeneratedNodeAndNewInputs, LinkedInputPort}
import wdl4s.wom.graph.GraphNodePort.{GraphNodeOutputPort, OutputPort}

sealed abstract class CallNode extends GraphNode {
  def callable: Callable
  def callType: String

  def portBasedInputs: Set[GraphNodePort.InputPort]
  def expressionBasedInputs: Map[String, InstantiatedExpression]

  override final val inputPorts = portBasedInputs ++ expressionBasedInputs.values.flatMap(_.inputPorts)
}

final case class TaskCallNode private(override val name: String, callable: TaskDefinition, portBasedInputs: Set[GraphNodePort.InputPort], expressionBasedInputs: Map[String, InstantiatedExpression]) extends CallNode {
  val callType: String = "task"
  override val outputPorts: Set[GraphNodePort.OutputPort] = {
    callable.outputs.map(o => GraphNodeOutputPort(o.name, o.womType, this))
  }
}

final case class WorkflowCallNode private(override val name: String, callable: WorkflowDefinition, portBasedInputs: Set[GraphNodePort.InputPort], expressionBasedInputs: Map[String, InstantiatedExpression]) extends CallNode {
  val callType: String = "workflow"
  override val outputPorts: Set[GraphNodePort.OutputPort] = {
    callable.innerGraph.nodes.collect { case gon: GraphOutputNode => GraphNodeOutputPort(gon.name, gon.womType, this) }
  }
}

object TaskCall {
  def graphFromDefinition(taskDefinition: TaskDefinition): ErrorOr[Graph] = {

    def linkOutput(call: GraphNode)(output: OutputDefinition): ErrorOr[GraphNode] = call.outputByName(output.name).map(out => PortBasedGraphOutputNode(output.name, output.womType, out))

    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap

    for {
      callWithInputs <- CallNode.callWithInputs(taskDefinition.name, taskDefinition, Map.empty, Set.empty, prefixSeparator = taskDefinition.prefixSeparator)
      outputs <- taskDefinition.outputs.toList.traverse(linkOutput(callWithInputs.node) _)
      callSet = Set[GraphNode](callWithInputs.node)
      inputsSet = callWithInputs.newInputs.toSet[GraphNode]
      outputsSet = outputs.toSet[GraphNode]
      graph <- Graph.validateAndConstruct(callSet ++ inputsSet ++ outputsSet)
    } yield graph
  }
}

object CallNode {

  final case class CallNodeAndNewInputs(node: CallNode, newInputs: Set[GraphInputNode]) extends GeneratedNodeAndNewInputs {
    def nodes: Set[GraphNode] = Set(node) ++ newInputs
  }

  /**
    * Don't use this directly; go via callWithInputs to make sure everything's in order when constructing a CallNode.
    */
  private[graph] def apply(name: String, callable: Callable, linkedInputPorts: Set[GraphNodePort.InputPort], instantiatedExpressionInputs: Map[String, InstantiatedExpression]): CallNode = callable match {
    case t: TaskDefinition => TaskCallNode(name, t, linkedInputPorts, instantiatedExpressionInputs)
    case w: WorkflowDefinition => WorkflowCallNode(name, w, linkedInputPorts, instantiatedExpressionInputs)
  }

  /**
    * Create a CallNode for a Callable (task or workflow).
    *
    * If an input is supplied as a port from another Node, it gets wired in directly.
    * If an input is supplied as an expression, we try to create an InstantiatedExpression and include that in the call.
    * If an input is not supplied, it gets created as a GraphInputNode.
    *
    */
  def callWithInputs(name: String, callable: Callable,
                     portInputs: Map[String, OutputPort],
                     expressionInputs: Set[GraphNodeInputExpression],
                     prefixSeparator: String = "."): ErrorOr[CallNodeAndNewInputs] = {

    val graphNodeSetter = new GraphNode.GraphNodeSetter()

    val instantiatedExpressionInputsAttempt: ErrorOr[Map[String, InstantiatedExpression]] = expressionInputs.toList traverse { _.instantiateExpressionWithInputName(graphNodeSetter) } map { _.toMap }

    instantiatedExpressionInputsAttempt map { instantiatedExpressionInputs =>
      val inputPortLinker = GraphNode.linkInputPort(name + prefixSeparator, portInputs, graphNodeSetter.get) _

      // Filter out the inputs we already have from expressions:
      val asYetUnsuppliedInputs = callable.inputs.filterNot(inputDef => instantiatedExpressionInputs.contains(inputDef.name))
      val linkedInputPortsAndGraphInputNodes = asYetUnsuppliedInputs map inputPortLinker

      val linkedInputPorts = linkedInputPortsAndGraphInputNodes.map(_.newInputPort)
      val graphInputNodes = linkedInputPortsAndGraphInputNodes collect { case LinkedInputPort(_, Some(gin)) => gin }

      val callNode = CallNode(name, callable, linkedInputPorts.toSet, instantiatedExpressionInputs)

      graphNodeSetter._graphNode = callNode
      CallNodeAndNewInputs(callNode, graphInputNodes.toSet)
    }
  }
}
