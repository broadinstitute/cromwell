package wdl.transforms.base.ast2wdlom

import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.{DeclarationContent, ExpressionElement, TypeElement}
import cats.syntax.apply._
import cats.syntax.either._
import common.transforms.CheckedAtoB

object AstToDeclarationContent {
  def astToDeclarationContent(implicit astNodeToTypeElement: CheckedAtoB[GenericAstNode, TypeElement],
                              astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement]
                             ): CheckedAtoB[GenericAst, DeclarationContent] = CheckedAtoB.fromErrorOr("read declaration") { a =>


    val nameValidation: ErrorOr[String] = astNodeToString(a.getAttribute("name")).toValidated
    val outputTypeValidation: ErrorOr[TypeElement] = a.getAttributeAs[TypeElement]("type").toValidated
    val expressionElementValidation: ErrorOr[ExpressionElement] = a.getAttributeAs[ExpressionElement]("expression").toValidated

    (nameValidation, outputTypeValidation, expressionElementValidation) mapN {
      (name, outputType, expression) => DeclarationContent(outputType, name, expression)
    }
  }
}
