package wom.graph

import lenthall.validation.ErrorOr.ErrorOr
import wdl.types.WdlType
import wom.expression.WomExpression
import wom.graph.GraphNodePort.{ConnectedInputPort, OutputPort}

sealed trait GraphOutputNode extends GraphNode {
  def womType: WdlType

  final override val outputPorts: Set[GraphNodePort.OutputPort] = Set.empty
}

/**
  * Exposes an existing output port as a graph output.
  */
final case class PortBasedGraphOutputNode(override val name: String, womType: WdlType, source: OutputPort) extends GraphOutputNode {
  val singleInputPort = ConnectedInputPort(name, womType, source, _ => this)
  override val inputPorts: Set[GraphNodePort.InputPort] = Set(singleInputPort)
}

/**
  * A graph output which is produced by evaluating an expression.
  *
  * NB: Construct this via ExpressionBasedGraphOutputNode.linkWithInputs(...)
  */
final case class ExpressionBasedGraphOutputNode private(override val name: String, instantiatedExpression: InstantiatedExpression) extends GraphOutputNode {
  override val womType = instantiatedExpression.womReturnType
  override val inputPorts = instantiatedExpression.inputPorts
}

object ExpressionBasedGraphOutputNode {
  def linkWithInputs(name: String, expression: WomExpression, inputMapping: Map[String, OutputPort]): ErrorOr[ExpressionBasedGraphOutputNode] =
    InstantiatedExpression.instantiateExpressionForNode(ExpressionBasedGraphOutputNode.apply)(name, expression, inputMapping)
}

