package wdl.transforms.base.ast2wdlom

import cats.syntax.apply._
import cats.syntax.validated._
import cats.syntax.either._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.{CommandPartElement, ExpressionElement, PlaceholderAttributeSet}
import wdl.model.draft3.elements.CommandPartElement.{PlaceholderCommandPartElement, StringCommandPartElement}

object AstNodeToCommandPartElement {
  def astNodeToCommandPartElement(implicit
    astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement],
    astNodeToPlaceholderAttributeSet: CheckedAtoB[GenericAstNode, PlaceholderAttributeSet]
  ): CheckedAtoB[GenericAstNode, CommandPartElement] = CheckedAtoB.fromErrorOr { a: GenericAstNode =>
    a match {
      case t: GenericTerminal => astNodeToString(t).toValidated map StringCommandPartElement
      case a: GenericAst =>
        val expressionElementV: ErrorOr[ExpressionElement] = a.getAttributeAs[ExpressionElement]("expr").toValidated
        val attributesV: ErrorOr[PlaceholderAttributeSet] =
          a.getAttributeAs[PlaceholderAttributeSet]("attributes").toValidated

        (expressionElementV, attributesV) mapN { (expressionElement, attributes) =>
          PlaceholderCommandPartElement(expressionElement, attributes)
        }
      case other => s"Conversion for $other not supported".invalidNel
    }
  }
}
