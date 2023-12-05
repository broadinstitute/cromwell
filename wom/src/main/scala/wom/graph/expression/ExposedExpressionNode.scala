package wom.graph.expression

import common.validation.ErrorOr.ErrorOr
import wom.expression.WomExpression
import wom.graph.GraphNodePort.{InputPort, OutputPort}
import wom.graph.WomIdentifier
import wom.types.WomType

object ExposedExpressionNode {
  def fromInputMapping(identifier: WomIdentifier,
                       expression: WomExpression,
                       explicitWomType: WomType,
                       inputMapping: Map[String, OutputPort]
  ): ErrorOr[ExposedExpressionNode] = {
    // This constructor ignores the evaluated type and uses the explicit type instead
    def constructor(identifier: WomIdentifier,
                    expression: WomExpression,
                    evaluatedType: WomType,
                    inputPorts: Map[String, InputPort]
    ) =
      new ExpressionNode(identifier, expression, explicitWomType, inputPorts) with ExposedExpressionNode
    ExpressionNode.buildFromConstructor(constructor)(identifier, expression, inputMapping)
  }
}

/**
  * An expression node that has an explicit WomType, that could be different from the evaluated type of the expression.
  * Coercion to this explicit type should be applied when evaluating the expression
  */
trait ExposedExpressionNode extends ExpressionNode
