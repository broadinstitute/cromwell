package wom.graph.expression

import lenthall.validation.ErrorOr.ErrorOr
import wom.expression.WomExpression
import wom.graph.GraphNodePort.{InputPort, OutputPort}
import wom.graph.WomIdentifier
import wom.types.WdlType

object AnonymousExpressionNode {
  def fromInputMapping(identifier: WomIdentifier,
                       expression: WomExpression,
                       inputMapping: Map[String, OutputPort]): ErrorOr[ExpressionNode] = {
    def constructor(identifier: WomIdentifier,
                    expression: WomExpression,
                    evaluatedType: WdlType,
                    inputPorts: Map[String, InputPort]) = {
      new ExpressionNode(identifier, expression, evaluatedType, inputPorts) with AnonymousExpressionNode
    }
    ExpressionNode.buildFromConstructor(constructor)(identifier, expression, inputMapping)
  }
}

/**
  * An expression node that is purely an internal expression and shouldn't be visible outside the graph
  */
trait AnonymousExpressionNode extends ExpressionNode
