package wdl4s.wom.graph

import cats.data.Validated.Valid
import cats.syntax.validated._
import cats.syntax.apply._
import lenthall.validation.ErrorOr.ErrorOr
import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
import wdl4s.wom.graph.GraphNodePort.{ConditionalOutputPort, InputPort, OutputPort}
import wdl4s.wdl.types.{WdlBooleanType, WdlOptionalType}
import wdl4s.wom.graph.GraphNode.{GeneratedNodeAndNewInputs, LinkedInputPort}

/**
  * Currently only WDL has the concept of conditional executions:
  *
  * @param innerGraph Imagine that the contents of the conditional block were a self-contained workflow. That's this Graph
  * @param condition The (boolean) expression on which the conditional is predicated.
  * @param innerGraphInputs The set of input ports that will be needed to satisfy the inner graph's inputs.
  * @param outputMapping Output ports for the conditional node which link back to GraphOutputNodes of the inner graph.
  */
final case class ConditionalNode private(innerGraph: Graph,
                                         condition: InstantiatedExpression,
                                         innerGraphInputs: Set[InputPort],
                                         outputMapping: Set[ConditionalOutputPort]) extends GraphNode {

  override val name: String = "ConditionalNode"

  override val inputPorts: Set[InputPort] = condition.inputPorts ++ innerGraphInputs
  override val outputPorts: Set[GraphNodePort.OutputPort] = outputMapping.toSet[OutputPort]
}

object ConditionalNode  {

  final case class ConditionalNodeWithInputs(node: ConditionalNode, newInputs: Set[GraphInputNode]) extends GeneratedNodeAndNewInputs {
    override val newExpressions = Set.empty[ExpressionNode]
  }

  def wireInConditional(innerGraph: Graph, condition: GraphNodeInputExpression, inputMapping: Map[String, OutputPort]): ErrorOr[ConditionalNodeWithInputs] = {
    val graphNodeSetter = new GraphNode.GraphNodeSetter()
    val inputPortLinker = GraphNode.linkInputPort("", inputMapping, graphNodeSetter.get) _

    val conditionalExpressionTypeValidation: ErrorOr[Unit] = condition.evaluateType flatMap {
      case WdlBooleanType => Valid(())
      case other => s"Cannot base a conditional on non-boolean type ${other.toWdlString}".invalidNel
    }

    // Try to instantiate the expression (i.e. connect its inputs with previous output ports)
    val conditionalExpressionValidation: ErrorOr[InstantiatedExpression] = condition.instantiateExpression(graphNodeSetter)

    // Create an input port for every input going into the ConditionalNode. If that meant creating a new GraphInputNode, that's
    // part of the LinkedInputPort case class.
    val linkedInputPortsAndGraphInputNodes: Set[LinkedInputPort] = innerGraph.nodes.inputDefinitions.map(inputPortLinker)
    // The next two lines split the LinkedInputPort case classes into the sets of InputPorts and GraphInputNodes
    val linkedInputPorts: Set[InputPort] = linkedInputPortsAndGraphInputNodes.map(_.newInputPort)
    val graphInputNodes: Set[GraphInputNode] = linkedInputPortsAndGraphInputNodes collect { case LinkedInputPort(_, Some(gin)) => gin }

    val outputPorts: Set[ConditionalOutputPort] = innerGraph.nodes.collect { case gon: PortBasedGraphOutputNode =>
      ConditionalOutputPort(gon.name, WdlOptionalType(gon.womType), gon, graphNodeSetter.get)
    }

    (conditionalExpressionValidation, conditionalExpressionTypeValidation) mapN { (conditionalExpression, _) =>
      val conditionalNode: ConditionalNode = ConditionalNode(innerGraph, conditionalExpression, linkedInputPorts, outputPorts)
      graphNodeSetter._graphNode = conditionalNode

      ConditionalNodeWithInputs(conditionalNode, graphInputNodes)
    }
  }
}
