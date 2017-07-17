package wdl4s.wom.graph

import wdl4s.wom.graph.GraphNodePort.{ConnectedInputPort, InputPort, OutputPort, ScatterGathererPort}
import wdl4s.wom.graph.ScatterNode.ScatterVariableMapping
import wdl4s.wdl.types.WdlArrayType
import wdl4s.wom.graph.GraphNode.LinkedInputPort

/**
  *
  * @param graph Imagine that the contents of a WDL scatter block were a self-contained workflow. That's this Graph
  * @param scatterVariableMapping Case class representing:
  *                               - An (Array[X] type) InputPort to be exposed by the ScatterNode
  *                               - The (X type) GraphInputNode in the inner graph which it provides the input for
  * @param otherInputPorts - Imagine the scatter were to be run as a workflow and you were writing an input JSON, these
  *                        are the inputs you might provide.
  *                        - NB It's 'other' because it doesn't include the scatter input variable.
  * @param outputMapping Output ports for the scatter node, which also link back to GraphOutputNodes of the inner graph.
  */
final case class ScatterNode private(graph: Graph,
                                     scatterVariableMapping: ScatterVariableMapping,
                                     otherInputPorts: Set[InputPort],
                                     outputMapping: Set[ScatterGathererPort]) extends GraphNode {

  def innerGraphInputs: Set[GraphInputNode] = graph.nodes.collect { case gin: GraphInputNode => gin }

  // NB if you find yourself calling .filter on this set of inputPorts, you probably just wanted to access either
  // the scatterVariableMapping or otherInputPorts fields directly.
  override val inputPorts: Set[InputPort] = Set(scatterVariableMapping.inputPort) ++ otherInputPorts
  override val outputPorts: Set[GraphNodePort.OutputPort] = outputMapping.toSet[OutputPort]
}

object ScatterNode {

  /**
    * Maps from an InputPort on the ScatterNode to the GraphInputNode in the innerGraph
    */
  case class ScatterVariableMapping(inputPort: InputPort, graphInputNode: GraphInputNode)

  final case class ScatterNodeWithInputs(scatter: ScatterNode, inputs: Set[GraphInputNode])

  /**
    * Creates a ScatterNode that's been connected to the inputs provided.
    * @param innerGraph The inner graph that's being scattered over.
    * @param scatterVariableSource The OutputPort (in the outer Graph) to connect to the ScatterNode's scatterVariable InputPort.
    * @param scatterVariableInput the GraphInputNode in the innerGraph which represents the scatter variable.
    * @param inputMapping The mapping to inputs expected by the innerGraph
    * @return The ScatterNodeWithInputs containing:
    *         * scatter: The scatter node itself, connected to the inputs
    *         * inputs: A set of GraphInputNodes for inputs which weren't supplied
    */
  // TODO: Validate that the scatter source is an array type (oh how I wish ports could be typed!)
  def scatterOverGraph(innerGraph: Graph, scatterVariableSource: OutputPort, scatterVariableInput: GraphInputNode, inputMapping: Map[String, OutputPort]): ScatterNodeWithInputs = {
    val graphNodeSetter = new GraphNode.GraphNodeSetter()
    val inputPortLinker = GraphNode.linkInputPort("", inputMapping, graphNodeSetter.get) _

    val scatterVariableMapping = ScatterVariableMapping(ConnectedInputPort(scatterVariableInput.name, scatterVariableInput.womType, scatterVariableSource, graphNodeSetter.get), scatterVariableInput)

    // TODO: Don't want to re-link the scatter variable...
    val linkedInputPortsAndGraphInputNodes: Set[LinkedInputPort] = (innerGraph.nodes - scatterVariableInput).inputDefinitions.map(inputPortLinker)
    val linkedInputPorts = linkedInputPortsAndGraphInputNodes.map(_.newInputPort)
    val graphInputNodes = linkedInputPortsAndGraphInputNodes collect { case LinkedInputPort(_, Some(gin)) => gin }

    val outputPorts: Set[ScatterGathererPort] = innerGraph.nodes.collect { case gon: GraphOutputNode =>
      ScatterGathererPort(gon.name, WdlArrayType(gon.womType), gon, graphNodeSetter.get)
    }

    val scatterNode: ScatterNode = ScatterNode(innerGraph, scatterVariableMapping, linkedInputPorts, outputPorts)
    graphNodeSetter._graphNode = scatterNode

    ScatterNodeWithInputs(scatterNode, graphInputNodes)
  }
}
