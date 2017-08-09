package wdl4s.wom.graph

import cats.implicits._
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.wom.callable.Callable.OutputDefinition
import wdl4s.wom.callable.{Callable, TaskDefinition, WorkflowDefinition}
import wdl4s.wom.graph.CallNode.CallWithInputs
import wdl4s.wom.graph.GraphNode.LinkedInputPort
import wdl4s.wom.graph.GraphNodePort.{GraphNodeOutputPort, OutputPort}

sealed abstract class CallNode extends GraphNode {
  def name: String
  def callable: Callable
  def callType: String
}

final case class TaskCallNode private(name: String, callable: TaskDefinition, inputPorts: Set[GraphNodePort.InputPort]) extends CallNode {
  val callType: String = "task"
  override val outputPorts: Set[GraphNodePort.OutputPort] = {
    callable.outputs.map(o => GraphNodeOutputPort(o.name, o.womType, this))
  }
}

final case class WorkflowCallNode private(name: String, callable: WorkflowDefinition, inputPorts: Set[GraphNodePort.InputPort]) extends CallNode {
  val callType: String = "workflow"
  override val outputPorts: Set[GraphNodePort.OutputPort] = {
    callable.innerGraph.nodes.collect { case gon: PortBasedGraphOutputNode => GraphNodeOutputPort(gon.name, gon.womType, this) }
  }
}

object TaskCall {
  def graphFromDefinition(taskDefinition: TaskDefinition): ErrorOr[Graph] = {

    def linkOutput(call: GraphNode)(output: OutputDefinition): ErrorOr[GraphNode] = call.outputByName(output.name).map(out => PortBasedGraphOutputNode(output.name, output.womType, out))

    val CallWithInputs(call, inputs) = CallNode.callWithInputs(taskDefinition.name, taskDefinition, Map.empty)
    val outputsValidation = taskDefinition.outputs.toList.traverse(linkOutput(call) _)

    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
    outputsValidation flatMap { outputs =>
      val callSet: Set[GraphNode] = Set[GraphNode](call)
      val inputsSet: Set[_ <: GraphNode] = inputs
      val outputsSet: Set[GraphNode] = outputs.toSet

      Graph.validateAndConstruct(callSet ++ inputsSet ++ outputsSet)
    }
  }
}

object CallNode {

  final case class CallWithInputs(call: CallNode, inputs: Set[GraphInputNode])

  /**
    * Create a CallNode for a Callable (task or workflow).
    *
    * If an input is supplied, it gets wired in as appropriate.
    * If an input is not supplied, it gets created as a GraphInputNode.
    *
    * The returned value is a tuple of (
    *   _1: the CallNode
    *   _2: any GraphInputNodes we created for unsupplied inputs
    * )
    */
  def callWithInputs(name: String, callable: Callable, inputMapping: Map[String, OutputPort]): CallWithInputs = {

    val graphNodeSetter = new GraphNode.GraphNodeSetter()
    val inputPortLinker = GraphNode.linkInputPort(callable.name + ".", inputMapping, graphNodeSetter.get) _
    val linkedInputPortsAndGraphInputNodes = callable.inputs map inputPortLinker

    val linkedInputPorts = linkedInputPortsAndGraphInputNodes.map(_.newInputPort)
    val graphInputNodes = linkedInputPortsAndGraphInputNodes collect { case LinkedInputPort(_, Some(gin)) => gin }

    val callNode = callable match {
      case t: TaskDefinition => TaskCallNode(name, t, linkedInputPorts)
      case w: WorkflowDefinition => WorkflowCallNode(name, w, linkedInputPorts)
    }

    graphNodeSetter._graphNode = callNode
    CallWithInputs(callNode, graphInputNodes)
  }
}
