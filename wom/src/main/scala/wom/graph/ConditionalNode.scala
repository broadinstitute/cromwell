package wom.graph

import wom.graph.GraphNode.GeneratedNodeAndNewNodes
import wom.graph.GraphNodePort.{ConditionalOutputPort, InputPort, OutputPort}
import wom.types.WdlOptionalType

/**
  * Currently only WDL has the concept of conditional executions:
  *
  * @param innerGraph Imagine that the contents of the conditional block were a self-contained workflow. That's this Graph
  * @param condition The (boolean) expression on which the conditional is predicated.
  * @param conditionalOutputPorts Output ports for the conditional node which link back to GraphOutputNodes of the inner graph.
  */
final case class ConditionalNode private(innerGraph: Graph,
                                         condition: ExpressionNode,
                                         conditionalOutputPorts: Set[ConditionalOutputPort]) extends GraphNode {

  override val identifier: WomIdentifier = WomIdentifier("ConditionalNode")

  override val inputPorts: Set[InputPort] = condition.inputPorts
  override val outputPorts: Set[GraphNodePort.OutputPort] = conditionalOutputPorts.toSet[OutputPort]
}

object ConditionalNode  {

  final case class ConditionalNodeWithNewNodes(node: ConditionalNode) extends GeneratedNodeAndNewNodes {
    override val newInputs = node.innerGraph.externalInputNodes.toSet[GraphInputNode]
    override val newExpressions = Set(node.condition)
  }

  def wireInConditional(innerGraph: Graph, expressionNode: ExpressionNode): ConditionalNodeWithNewNodes = {
    val graphNodeSetter = new GraphNode.GraphNodeSetter()

    val outputPorts: Set[ConditionalOutputPort] = innerGraph.nodes.collect { case gon: PortBasedGraphOutputNode =>
      ConditionalOutputPort(WdlOptionalType(gon.womType), gon, graphNodeSetter.get)
    }

    val conditionalNode: ConditionalNode = ConditionalNode(innerGraph, expressionNode, outputPorts)
    graphNodeSetter._graphNode = conditionalNode

    ConditionalNodeWithNewNodes(conditionalNode)
  }
}
