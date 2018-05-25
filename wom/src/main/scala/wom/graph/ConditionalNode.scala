package wom.graph

import wom.graph.GraphNode.{GeneratedNodeAndNewNodes, GraphNodeWithInnerGraph}
import wom.graph.GraphNodePort.{ConditionalOutputPort, ConnectedInputPort, InputPort, OutputPort}
import wom.graph.expression.ExpressionNode
import wom.types.WomBooleanType
import common.collections.EnhancedCollections._

/**
  * Currently only WDL has the concept of conditional executions:
  *
  * @param innerGraph Imagine that the contents of the conditional block were a self-contained workflow. That's this Graph
  * @param conditionExpression The (boolean) expression on which the conditional is predicated.
  * @param conditionalOutputPorts Output ports for the conditional node which link back to GraphOutputNodes of the inner graph.
  */
final case class ConditionalNode private(override val innerGraph: Graph,
                                         conditionExpression: ExpressionNode,
                                         conditionalOutputPorts: Set[ConditionalOutputPort]) extends GraphNode with GraphNodeWithInnerGraph {

  override val identifier: WomIdentifier = WomIdentifier("ConditionalNode")

  override val inputPorts: Set[InputPort] = Set(ConnectedInputPort("condition", WomBooleanType, conditionExpression.singleOutputPort, _ => this))
  override val outputPorts: Set[GraphNodePort.OutputPort] = conditionalOutputPorts.toSet[OutputPort]
}

object ConditionalNode  {

  final case class ConditionalNodeWithNewNodes(node: ConditionalNode) extends GeneratedNodeAndNewNodes {
    override val newInputs = node.innerGraph.externalInputNodes
    override val usedOuterGraphInputNodes =
      (node.conditionExpression.upstream.filterByType[OuterGraphInputNode]: Set[OuterGraphInputNode]) ++
        (node.innerGraph.outerGraphInputNodes.map(_.linkToOuterGraphNode).filterByType[OuterGraphInputNode]: Set[OuterGraphInputNode])

    override val newExpressions = Set(node.conditionExpression)
  }

  def wireInConditional(innerGraph: Graph, expressionNode: ExpressionNode): ConditionalNodeWithNewNodes = {
    val graphNodeSetter = new GraphNode.GraphNodeSetter[ConditionalNode]()

    val outputPorts: Set[ConditionalOutputPort] = innerGraph.nodes.collect { case gon: PortBasedGraphOutputNode =>
      ConditionalOutputPort(gon, graphNodeSetter.get)
    }

    val conditionalNode: ConditionalNode = ConditionalNode(innerGraph, expressionNode, outputPorts)
    graphNodeSetter._graphNode = conditionalNode

    ConditionalNodeWithNewNodes(conditionalNode)
  }
}
