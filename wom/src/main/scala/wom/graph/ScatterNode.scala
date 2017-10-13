package wom.graph

import wom.graph.GraphNode.GeneratedNodeAndNewNodes
import wom.graph.GraphNodePort.{InputPort, OutputPort, ScatterGathererPort}
import wom.types.WdlArrayType

/**
  *
  * @param innerGraph Imagine that the contents of a WDL scatter block were a self-contained workflow. That's this Graph
  * @param scatterCollectionExpressionNode expression node for the scatter collection expression
  * @param scatterVariableInnerGraphInputNode graph input node in the inner graph satisfied by the expression node
  * @param outputMapping Output ports for the scatter node, which also link back to GraphOutputNodes of the inner graph.
  */
final case class ScatterNode private(innerGraph: Graph,
                                     scatterCollectionExpressionNode: ExpressionNode,
                                     scatterVariableInnerGraphInputNode: GraphInputNode,
                                     outputMapping: Set[ScatterGathererPort]) extends GraphNode {

  override val identifier: WomIdentifier = WomIdentifier("ScatterNode")

  // NB if you find yourself calling .filter on this set of inputPorts, you probably just wanted to access either
  // the scatterVariableMapping or otherInputPorts fields directly.
  override val inputPorts: Set[InputPort] = scatterCollectionExpressionNode.inputPorts
  override val outputPorts: Set[GraphNodePort.OutputPort] = outputMapping.toSet[OutputPort]

  lazy val nodes: Set[GraphNode] = Set(this, scatterCollectionExpressionNode)
}

object ScatterNode {

  /**
    * Maps from an InstantiatedExpression on the ScatterNode to the GraphInputNode in the innerGraph.
    * A 'None' for graphInputNode indicates that the inner graph doesn't use the scatter variable.
    */
  case class ScatterVariableMapping(scatterCollectionExpressionNode: ExpressionNode, graphInputNode: GraphInputNode)

  final case class ScatterNodeWithNewNodes(node: ScatterNode) extends GeneratedNodeAndNewNodes {
    override val newExpressions = Set(node.scatterCollectionExpressionNode)
    override val newInputs = node.innerGraph.externalInputNodes.toSet[GraphInputNode]
  }

  /**
    * Helper class to build call nodes.
    * Helps making input ports and building the node while making sure node references are set properly.
    */
  def scatterOverGraph(innerGraph: Graph,
                       scatterCollectionExpressionNode: ExpressionNode,
                       scatterVariableInnerGraphInputNode: GraphInputNode): ScatterNodeWithNewNodes = {
    val graphNodeSetter = new GraphNode.GraphNodeSetter()

    val outputPorts: Set[ScatterGathererPort] = innerGraph.nodes.collect { case gon: PortBasedGraphOutputNode =>
      ScatterGathererPort(WdlArrayType(gon.womType), gon, graphNodeSetter.get)
    }

    val scatterNode = ScatterNode(
      innerGraph,
      scatterCollectionExpressionNode,
      scatterVariableInnerGraphInputNode,
      outputPorts
    )

    graphNodeSetter._graphNode = scatterNode

    ScatterNodeWithNewNodes(scatterNode)
  }
}
