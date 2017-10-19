package wom.graph.expression

import lenthall.validation.ErrorOr.ErrorOr
import wom.expression.WomExpression
import wom.graph.GraphNodePort.{InputPort, OutputPort}
import wom.graph.WomIdentifier
import wom.types.WdlType

object ExposedExpressionNode {
  def fromInputMapping(identifier: WomIdentifier,
                       expression: WomExpression,
                       explicitWomType: WdlType, 
                       inputMapping: Map[String, OutputPort]): ErrorOr[ExpressionNode] = {
    // This constructor ignores the evaluated type and uses the explicit type instead
    def constructor(identifier: WomIdentifier,
                    expression: WomExpression,
                    evaluatedType: WdlType,
                    inputPorts: Map[String, InputPort]) = {
      new ExpressionNode(identifier, expression, explicitWomType, inputPorts) with AnonymousExpressionNode
    }
    ExpressionNode.buildFromConstructor(constructor)(identifier, expression, inputMapping)
  }
}

/**
  * An expression node that has an explicit WomType, that could be different from the evaluated type of the expression.
  * Coercion to this explicit type should be applied when evaluating the expression
  */
trait ExposedExpressionNode extends ExpressionNode
