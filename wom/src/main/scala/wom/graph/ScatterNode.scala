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
  * @param scatterVariableNodes Inner graph nodes for each scatter collection expression being scattered over. WDL uses exactly one, CWL >= 1.
  * @param outputMapping Output ports for the scatter node, which also link back to GraphOutputNodes of the inner graph.
  */
final case class ScatterNode private (override val innerGraph: Graph,
                                      scatterVariableNodes: List[ScatterVariableNode],
                                      outputMapping: Set[ScatterGathererPort],
                                      scatterProcessingFunction: ScatterProcessingFunction,
                                      scatterCollectionFunctionBuilder: ScatterCollectionFunctionBuilder
) extends GraphNode
    with GraphNodeWithInnerGraph {

  override val identifier: WomIdentifier = WomIdentifier("ScatterNode")

  val scatterCollectionExpressionNodes: List[ExpressionNode] = scatterVariableNodes.map(_.scatterExpressionNode)
  val scatterVariableInnerGraphInputNodes: List[OuterGraphInputNode] = scatterVariableNodes

  // NB if you find yourself calling .filter on this set of inputPorts, you probably just wanted to access either
  // the scatterVariableMapping or otherInputPorts fields directly.
  override val inputPorts: Set[InputPort] = scatterCollectionExpressionNodes.toSet[ExpressionNode] map {
    scatterCollectionExpressionNode =>
      ConnectedInputPort(
        scatterCollectionExpressionNode.identifier.localName.value,
        scatterCollectionExpressionNode.womType,
        scatterCollectionExpressionNode.singleOutputPort,
        _ => this
      )
  }
  override val outputPorts: Set[GraphNodePort.OutputPort] = outputMapping.toSet[OutputPort]

  lazy val nodes: Set[GraphNode] = scatterCollectionExpressionNodes.toSet + this
}

object ScatterNode {
  case class ScatterVariableAndValue(scatterVariableNode: ScatterVariableNode, arrayValue: WomArray)

  type ScatterProcessingFunction = List[ScatterVariableAndValue] => Checked[Int]
  // Use dot product as the default processing function (works well for single variable scatters too)
  val DefaultScatterProcessingFunction: ScatterProcessingFunction = { nodesAndValues: List[ScatterVariableAndValue] =>
    nodesAndValues.map(_.arrayValue.size).distinct match {
      case head :: Nil => head.validNelCheck
      case _ =>
        "All arrays must have the same number of element when using the dot product scatter method".invalidNelCheck
    }
  }

  type ScatterCollectionFunction = (List[WomValue], WomArrayType) => WomArray
  val DefaultScatterCollectionFunction: ScatterCollectionFunction = {
    (shards: List[WomValue], valueType: WomArrayType) => WomArray(valueType, shards)
  }

  type ScatterCollectionFunctionBuilder = List[Int] => ScatterCollectionFunction
  val DefaultScatterCollectionFunctionBuilder: ScatterCollectionFunctionBuilder = { _: List[Int] =>
    DefaultScatterCollectionFunction
  }

  /**
    * Maps from an InstantiatedExpression on the ScatterNode to the GraphInputNode in the innerGraph.
    * A 'None' for graphInputNode indicates that the inner graph doesn't use the scatter variable.
    */
  case class ScatterVariableMapping(scatterCollectionExpressionNode: ExpressionNode, graphInputNode: GraphInputNode)

  final case class ScatterNodeWithNewNodes(node: ScatterNode) extends GeneratedNodeAndNewNodes {
    override val newExpressions: Set[ExpressionNode] = node.scatterCollectionExpressionNodes.toSet
    override val newInputs: Set[ExternalGraphInputNode] = node.innerGraph.externalInputNodes
    override val usedOuterGraphInputNodes: Set[OuterGraphInputNode] =
      (node.scatterCollectionExpressionNodes
        .flatMap(_.upstream)
        .toSet
        .filterByType[OuterGraphInputNode]: Set[OuterGraphInputNode]) ++
        (node.innerGraph.outerGraphInputNodes
          .map(_.linkToOuterGraphNode)
          .filterByType[OuterGraphInputNode]: Set[OuterGraphInputNode])
  }

  /**
    * Helper class to build call nodes.
    * Helps making input ports and building the node while making sure node references are set properly.
    */
  def scatterOverGraph(innerGraph: Graph,
                       scatterVariableInnerGraphInputNode: ScatterVariableNode
  ): ScatterNodeWithNewNodes = {
    val scatterNodeBuilder = new ScatterNodeBuilder
    val outputPorts: Set[ScatterGathererPort] = innerGraph.nodes.collect { case gon: PortBasedGraphOutputNode =>
      scatterNodeBuilder.makeOutputPort(WomArrayType(gon.womType), gon)
    }

    scatterNodeBuilder.build(innerGraph,
                             outputPorts,
                             List(scatterVariableInnerGraphInputNode),
                             DefaultScatterProcessingFunction,
                             DefaultScatterCollectionFunctionBuilder
    )
  }

  /**
    * Helper class to build call nodes.
    * Helps making output ports and building the node while making sure node references are set properly.
    */
  class ScatterNodeBuilder {
    private val graphNodeSetter = new GraphNode.GraphNodeSetter[ScatterNode]()

    def makeOutputPort(womType: WomArrayType, nodeToGather: PortBasedGraphOutputNode): ScatterGathererPort =
      ScatterGathererPort(womType, nodeToGather, graphNodeSetter.get)

    def build(innerGraph: Graph,
              outputPorts: Set[ScatterGathererPort],
              scatterVariableNodes: List[ScatterVariableNode],
              scatterProcessingFunction: ScatterProcessingFunction,
              scatterCollectionFunctionBuilder: ScatterCollectionFunctionBuilder
    ) = {
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
}
