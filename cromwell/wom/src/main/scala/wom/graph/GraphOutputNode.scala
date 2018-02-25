package wom.graph

import common.validation.ErrorOr.ErrorOr
import wom.expression.WomExpression
import wom.graph.GraphNodePort.{ConnectedInputPort, InputPort, OutputPort}
import wom.graph.expression.ExpressionNode
import wom.graph.expression.ExpressionNode.buildFromConstructor
import wom.types.WomType

sealed trait GraphOutputNode extends GraphNode {
  def womType: WomType
  def graphOutputPort: OutputPort
}

/**
  * Exposes an existing output port as a graph output.
  */
final case class PortBasedGraphOutputNode(override val identifier: WomIdentifier, womType: WomType, source: OutputPort) extends GraphOutputNode {
  val singleInputPort: InputPort = ConnectedInputPort(localName, womType, source, _ => this)
  lazy val singleUpstreamNode: GraphNode = singleInputPort.upstream.graphNode
  lazy val singleUpstreamPort: OutputPort = singleInputPort.upstream
  override val inputPorts: Set[GraphNodePort.InputPort] = Set(singleInputPort)
  override val outputPorts: Set[GraphNodePort.OutputPort] = Set.empty
  override val graphOutputPort: OutputPort = singleInputPort.upstream
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
trait ExpressionBasedGraphOutputNode extends ExpressionNode with GraphOutputNode {
  override val graphOutputPort: OutputPort = singleOutputPort
}
