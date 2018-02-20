package wdl.draft3.transforms.ast2wdlom

import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser
import wdl.model.draft3.elements.{DeclarationContent, ExpressionElement, TypeElement}
import cats.syntax.apply._
import cats.syntax.either._
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._

object AstToDeclarationContent {
  def convert(a: WdlParser.Ast): ErrorOr[DeclarationContent] = {

    val nameValidation: ErrorOr[String] = astNodeToString(a.getAttribute("name")).toValidated
    val outputTypeValidation: ErrorOr[TypeElement] = a.getAttributeAs[TypeElement]("type").toValidated
    val expressionElementValidation: ErrorOr[ExpressionElement] = a.getAttributeAs[ExpressionElement]("expression").toValidated

    (nameValidation, outputTypeValidation, expressionElementValidation) mapN {
      (name, outputType, expression) => DeclarationContent(outputType, name, expression)
    }
  }
}
