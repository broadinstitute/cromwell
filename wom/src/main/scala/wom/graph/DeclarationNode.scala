package wom.graph

import lenthall.validation.ErrorOr.ErrorOr
import wdl.types.WdlType
import wom.expression.WomExpression
import wom.graph.ExpressionNode.buildFromConstructor
import wom.graph.GraphNodePort.{InputPort, OutputPort}

object DeclarationNode {
  def fromInputMapping(identifier: WomIdentifier,
                       expression: WomExpression,
                       explicitWomType: WdlType, 
                       inputMapping: Map[String, OutputPort]): ErrorOr[ExpressionNode] = {
    // This constructor ignores the evaluated type and uses the explicit type instead
    def constructor(identifier: WomIdentifier,
                    expression: WomExpression,
                    evaluatedType: WdlType,
                    inputPorts: Map[String, InputPort]) = {
      new ExpressionNode(identifier, expression, explicitWomType, inputPorts) with DeclarationNode
    }
    buildFromConstructor(constructor)(identifier, expression, inputMapping)
  }
}

/**
  * An expression node that has an explicit WomType, that could be different from the evaluated type of the expression.
  * Coercion to this explicit type should be applied when evaluating the expression
  */
trait DeclarationNode extends ExpressionNode
