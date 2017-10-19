package wom.graph

import lenthall.validation.ErrorOr.ErrorOr
import wdl.types.WdlType
import wom.expression.WomExpression
import wom.graph.expression.ExpressionNode.buildFromConstructor
import wom.graph.GraphNodePort.{ConnectedInputPort, InputPort, OutputPort}
import wom.graph.expression.{AnonymousExpressionNode, ExpressionNode}

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

object ExpressionBasedGraphOutputNode {
  def fromInputMapping(identifier: WomIdentifier,
                       expression: WomExpression,
                       explicitWomType: WdlType,
                       inputMapping: Map[String, OutputPort]): ErrorOr[ExpressionBasedGraphOutputNode] = {
    // This constructor ignores the evaluated type and uses the explicit type instead
    def constructor(identifier: WomIdentifier,
                    expression: WomExpression,
                    evaluatedType: WdlType,
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
