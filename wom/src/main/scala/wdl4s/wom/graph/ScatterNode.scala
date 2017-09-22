package wdl4s.wom.graph

import cats.data.Validated.Valid
import cats.syntax.validated._
import cats.syntax.apply._
import lenthall.validation.ErrorOr.ErrorOr
import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
import wdl4s.wom.graph.GraphNodePort.{InputPort, OutputPort, ScatterGathererPort}
import wdl4s.wom.graph.ScatterNode.ScatterVariableMapping
import wdl4s.wdl.types.{WdlArrayType, WdlMapType}
import wdl4s.wom.graph.GraphNode.{GeneratedNodeAndNewInputs, LinkedInputPort}

/**
  *
  * @param innerGraph Imagine that the contents of a WDL scatter block were a self-contained workflow. That's this Graph
  * @param scatterVariableMapping Case class representing:
  *                               - An (Array[X] type) InputPort to be exposed by the ScatterNode
  *                               - The (X type) GraphInputNode in the inner graph which it provides the input for
  * @param otherInputPorts - Imagine the scatter were to be run as a workflow and you were writing an input JSON, these
  *                        are the inputs you might provide.
  *                        - NB It's 'other' because it doesn't include the scatter input variable.
  * @param outputMapping Output ports for the scatter node, which also link back to GraphOutputNodes of the inner graph.
  */
final case class ScatterNode private(innerGraph: Graph,
                                     scatterVariableMapping: ScatterVariableMapping,
                                     otherInputPorts: Set[InputPort],
                                     outputMapping: Set[ScatterGathererPort]) extends GraphNode {

  override val name: String = "ScatterNode"

  // NB if you find yourself calling .filter on this set of inputPorts, you probably just wanted to access either
  // the scatterVariableMapping or otherInputPorts fields directly.
  override val inputPorts: Set[InputPort] = scatterVariableMapping.scatterInstantiatedExpression.inputPorts ++ otherInputPorts
  override val outputPorts: Set[GraphNodePort.OutputPort] = outputMapping.toSet[OutputPort]
}

object ScatterNode {

  /**
    * Maps from an InstantiatedExpression on the ScatterNode to the GraphInputNode in the innerGraph.
    * A 'None' for graphInputNode indicates that the inner graph doesn't use the scatter variable.
    */
  case class ScatterVariableMapping(scatterInstantiatedExpression: InstantiatedExpression, graphInputNode: GraphInputNode)

  final case class ScatterNodeWithInputs(node: ScatterNode, newInputs: Set[GraphInputNode]) extends GeneratedNodeAndNewInputs {
    override val newExpressions = Set.empty[ExpressionNode]
  }

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
  // TODO: For CWL, we'll probably want the 'scatterVariableSource' to be an expression OR a GraphNodePort.OutputPort. TBC...
  def scatterOverGraph(innerGraph: Graph, scatterVariableSource: GraphNodeInputExpression, scatterVariableInput: GraphInputNode, inputMapping: Map[String, OutputPort]): ErrorOr[ScatterNodeWithInputs] = {
    val graphNodeSetter = new GraphNode.GraphNodeSetter()
    val inputPortLinker = GraphNode.linkInputPort("", inputMapping, graphNodeSetter.get) _

    val scatterExpressionType: ErrorOr[WdlArrayType] = scatterVariableSource.evaluateType flatMap {
      case arrayType: WdlArrayType => Valid(arrayType)
      case mapType: WdlMapType => Valid(mapType.equivalentArrayType)
      case other => s"Cannot scatter over non-array type ${other.toWdlString}".invalidNel
    }

    val scatterVariableMappingValidation: ErrorOr[ScatterVariableMapping] = scatterVariableSource.instantiateExpression(graphNodeSetter) map { ScatterVariableMapping(_, scatterVariableInput) }

    // Filter because we don't want to re-link the scatter variable...
    val linkedInputPortsAndGraphInputNodes: Set[LinkedInputPort] = (innerGraph.nodes - scatterVariableInput).inputDefinitions.map(inputPortLinker)
    val linkedInputPorts = linkedInputPortsAndGraphInputNodes.map(_.newInputPort)
    val graphInputNodes = linkedInputPortsAndGraphInputNodes collect { case LinkedInputPort(_, Some(gin)) => gin }

    val outputPorts: Set[ScatterGathererPort] = innerGraph.nodes.collect { case gon: PortBasedGraphOutputNode =>
      ScatterGathererPort(gon.name, WdlArrayType(gon.womType), gon, graphNodeSetter.get)
    }

    (scatterVariableMappingValidation, scatterExpressionType) mapN { (scatterVariableMapping, _) =>
      val scatterNode: ScatterNode = ScatterNode(innerGraph, scatterVariableMapping, linkedInputPorts, outputPorts)
      graphNodeSetter._graphNode = scatterNode

      ScatterNodeWithInputs(scatterNode, graphInputNodes)
    }
  }
}
