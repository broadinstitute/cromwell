package wdl.draft3.transforms.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.Ast
import wdl.model.draft3.elements.{ExpressionElement, InputDeclarationElement, TypeElement}

object AstToInputDeclarationElement {
  def convert(a: Ast): ErrorOr[InputDeclarationElement] = {

    val nameValidation: ErrorOr[String] = astNodeToString(a.getAttribute("name")).toValidated
    val inputTypeValidation: ErrorOr[TypeElement] = astNodeToTypeElement(a.getAttribute("type")).toValidated
    val expressionValidation: ErrorOr[Option[ExpressionElement]] = Option(a.getAttribute("expression")) match {
      case Some(e) => astNodeToExpressionElement(e).map(Some.apply).toValidated
      case None => None.validNel
    }

    (nameValidation, inputTypeValidation, expressionValidation) mapN { (name, inputType, expression) => InputDeclarationElement(inputType, name, expression) }
  }
}
