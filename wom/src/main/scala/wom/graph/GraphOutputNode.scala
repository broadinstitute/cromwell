package wom.graph

import lenthall.validation.ErrorOr.ErrorOr
import wom.expression.WomExpression
import wom.graph.GraphNodePort.{ConnectedInputPort, InputPort, OutputPort}
import wom.graph.expression.ExpressionNode.buildFromConstructor
import wom.graph.expression.{AnonymousExpressionNode, ExpressionNode}
import wom.types.WomType

sealed trait GraphOutputNode extends GraphNode {
  def womType: WomType
}

/**
  * Exposes an existing output port as a graph output.
  */
final case class PortBasedGraphOutputNode(override val identifier: WomIdentifier, womType: WomType, source: OutputPort) extends GraphOutputNode {
  val singleInputPort = ConnectedInputPort(localName, womType, source, _ => this)
  override val inputPorts: Set[GraphNodePort.InputPort] = Set(singleInputPort)
  override val outputPorts: Set[GraphNodePort.OutputPort] = Set(source)
}

object ExpressionBasedGraphOutputNode {
  def fromInputMapping(identifier: WomIdentifier,
                       expression: WomExpression,
                       explicitWomType: WomType,
                       inputMapping: Map[String, OutputPort]): ErrorOr[ExpressionBasedGraphOutputNode] = {
    // This constructor ignores the evaluated type and uses the explicit type instead
    def constructor(identifier: WomIdentifier,
                    expression: WomExpression,
                    evaluatedType: WomType,
                    inputPorts: Map[String, InputPort]) = {
      new ExpressionNode(identifier, expression, explicitWomType, inputPorts) with ExpressionBasedGraphOutputNode
    }
    buildFromConstructor(constructor)(identifier, expression, inputMapping)
  }
}

/**
  * A graph output which is produced by evaluating an expression.
  */
trait ExpressionBasedGraphOutputNode extends AnonymousExpressionNode with GraphOutputNode
