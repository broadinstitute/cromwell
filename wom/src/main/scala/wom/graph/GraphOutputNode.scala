package wom.graph

import lenthall.validation.ErrorOr.ErrorOr
import wom.expression.WomExpression
import wom.graph.GraphNodePort.{ConnectedInputPort, GraphNodeOutputPort, OutputPort}
import wom.types.WdlType

sealed trait GraphOutputNode extends GraphNode {
  def womType: WdlType
}

/**
  * Exposes an existing output port as a graph output.
  */
final case class PortBasedGraphOutputNode(override val identifier: WomIdentifier, womType: WdlType, source: OutputPort) extends GraphOutputNode {
  val singleInputPort = ConnectedInputPort(localName, womType, source, _ => this)
  override val inputPorts: Set[GraphNodePort.InputPort] = Set(singleInputPort)
  override val outputPorts: Set[GraphNodePort.OutputPort] = Set(source)
}

/**
  * A graph output which is produced by evaluating an expression.
  *
  * NB: Construct this via ExpressionBasedGraphOutputNode.linkWithInputs(...)
  */
final case class ExpressionBasedGraphOutputNode private(override val identifier: WomIdentifier, womType: WdlType, instantiatedExpression: InstantiatedExpression) extends GraphOutputNode {
  override val inputPorts = instantiatedExpression.inputPorts
  val singleOutputPort = GraphNodeOutputPort(identifier, womType, this)
  override val outputPorts: Set[GraphNodePort.OutputPort] = Set(singleOutputPort)
}

object ExpressionBasedGraphOutputNode {
  def linkWithInputs(nodeIdentifier: WomIdentifier, womType: WdlType, expression: WomExpression, inputMapping: Map[String, OutputPort]): ErrorOr[ExpressionBasedGraphOutputNode] = {
    val nodeConstructor = (identifier: WomIdentifier, instantiatedExpression: InstantiatedExpression) => {
      ExpressionBasedGraphOutputNode(identifier, womType, instantiatedExpression)
    }
    InstantiatedExpression.instantiateExpressionForNode(nodeConstructor)(nodeIdentifier, expression, inputMapping)
  }
}
