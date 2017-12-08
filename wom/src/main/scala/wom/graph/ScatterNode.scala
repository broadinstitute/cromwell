package wom.graph

import common.Checked
import common.collections.EnhancedCollections._
import common.validation.Checked._
import wom.graph.GraphNode.{GeneratedNodeAndNewNodes, GraphNodeWithInnerGraph}
import wom.graph.GraphNodePort.{ConnectedInputPort, InputPort, OutputPort, ScatterGathererPort}
import wom.graph.ScatterNode.{ScatterCollectionFunctionBuilder, ScatterProcessingFunction}
import wom.graph.expression.ExpressionNode
import wom.types.WomArrayType
import wom.values.{WomArray, WomValue}

/**
  *
  * @param innerGraph Imagine that the contents of a WDL scatter block were a self-contained workflow. That's this Graph
  * @param scatterVariableNodes inner graph nodes for each scatter collection expression being scattered over
  * @param outputMapping Output ports for the scatter node, which also link back to GraphOutputNodes of the inner graph.
  */
final case class ScatterNode private(override val innerGraph: Graph,
                                     scatterVariableNodes: List[ScatterVariableNode],
                                     outputMapping: Set[ScatterGathererPort],
                                     scatterProcessingFunction: ScatterProcessingFunction,
                                     scatterCollectionFunctionBuilder: ScatterCollectionFunctionBuilder) extends GraphNode with GraphNodeWithInnerGraph {

  override val identifier: WomIdentifier = WomIdentifier("ScatterNode")

  val scatterCollectionExpressionNodes: List[ExpressionNode] = scatterVariableNodes.map(_.scatterExpressionNode)
  val scatterVariableInnerGraphInputNodes: List[OuterGraphInputNode] = scatterVariableNodes

  // NB if you find yourself calling .filter on this set of inputPorts, you probably just wanted to access either
  // the scatterVariableMapping or otherInputPorts fields directly.
  override val inputPorts: Set[InputPort] = scatterCollectionExpressionNodes.toSet[ExpressionNode] map { scatterCollectionExpressionNode =>
    ConnectedInputPort(
      scatterCollectionExpressionNode.identifier.localName.value,
      scatterCollectionExpressionNode.womType,
      scatterCollectionExpressionNode.singleExpressionOutputPort,
      _ => this
    )
  }
  override val outputPorts: Set[GraphNodePort.OutputPort] = outputMapping.toSet[OutputPort]

  lazy val nodes: Set[GraphNode] = scatterCollectionExpressionNodes.toSet + this
}

object ScatterNode {
  case class ScatterVariableAndValue(scatterVariableNode: ScatterVariableNode, arrayValue: WomArray)
  
  type ScatterProcessingFunction = List[ScatterVariableAndValue] => Checked[Int]
  val DefaultScatterProcessingFunction: ScatterProcessingFunction = { scatterArrays: List[ScatterVariableAndValue] => 
    scatterArrays.map(_.arrayValue.size).sum.validNelCheck 
  }
  
  type ScatterCollectionFunction = (List[WomValue], WomArrayType) => WomArray
  val DefaultScatterCollectionFunction: ScatterCollectionFunction = { (shards: List[WomValue], valueType: WomArrayType) => WomArray(valueType, shards) }
  
  type ScatterCollectionFunctionBuilder = List[ScatterVariableAndValue] => ScatterCollectionFunction
  val DefaultScatterCollectionFunctionBuilder: ScatterCollectionFunctionBuilder = { _: List[ScatterVariableAndValue] => DefaultScatterCollectionFunction }

  /**
    * Maps from an InstantiatedExpression on the ScatterNode to the GraphInputNode in the innerGraph.
    * A 'None' for graphInputNode indicates that the inner graph doesn't use the scatter variable.
    */
  case class ScatterVariableMapping(scatterCollectionExpressionNode: ExpressionNode, graphInputNode: GraphInputNode)

  final case class ScatterNodeWithNewNodes(node: ScatterNode) extends GeneratedNodeAndNewNodes {
    override val newExpressions: Set[ExpressionNode] = node.scatterCollectionExpressionNodes.toSet
    override val newInputs: Set[ExternalGraphInputNode] = node.innerGraph.externalInputNodes
    override val usedOuterGraphInputNodes: Set[OuterGraphInputNode] =(node.scatterCollectionExpressionNodes.flatMap(_.upstream).toSet.filterByType[OuterGraphInputNode]: Set[OuterGraphInputNode]) ++
      (node.innerGraph.outerGraphInputNodes.map(_.linkToOuterGraphNode).filterByType[OuterGraphInputNode]: Set[OuterGraphInputNode])
    def nodes: Set[GraphNode] = newExpressions ++ newInputs ++ usedOuterGraphInputNodes ++ Set(node)
  }

  /**
    * Helper class to build call nodes.
    * Helps making input ports and building the node while making sure node references are set properly.
    */
  def scatterOverGraph(innerGraph: Graph,
                       scatterVariableInnerGraphInputNode: ScatterVariableNode): ScatterNodeWithNewNodes = {
    scatterOverGraph(innerGraph, List(scatterVariableInnerGraphInputNode), DefaultScatterProcessingFunction, DefaultScatterCollectionFunctionBuilder)
  }

  /**
    * Helper class to build call nodes.
    * Helps making input ports and building the node while making sure node references are set properly.
    */
  def scatterOverGraph(innerGraph: Graph,
                       scatterVariableNodes: List[ScatterVariableNode],
                       scatterProcessingFunction: ScatterProcessingFunction,
                       scatterCollectionFunctionBuilder: ScatterCollectionFunctionBuilder): ScatterNodeWithNewNodes = {
    val graphNodeSetter = new GraphNode.GraphNodeSetter[ScatterNode]()

    val outputPorts: Set[ScatterGathererPort] = innerGraph.nodes.collect { case gon: PortBasedGraphOutputNode =>
      ScatterGathererPort(WomArrayType(gon.womType), gon, graphNodeSetter.get)
    }

    val scatterNode = ScatterNode(
      innerGraph,
      scatterVariableNodes,
      outputPorts,
      scatterProcessingFunction,
      scatterCollectionFunctionBuilder
    )

    graphNodeSetter._graphNode = scatterNode

    ScatterNodeWithNewNodes(scatterNode)
  }
}
