package wdl.draft3.transforms.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.Ast
import wdl.model.draft3.elements.{InputDeclarationElement, TypeElement}

object AstToInputDeclarationElement {
  def convert(a: Ast): ErrorOr[InputDeclarationElement] = {

    val nameValidation: ErrorOr[String] = astNodeToString(a.getAttribute("name")).toValidated
    val inputTypeValidation: ErrorOr[TypeElement] = astNodeToTypeElement(a.getAttribute("type")).toValidated

    (nameValidation, inputTypeValidation) mapN { (name, inputType) => InputDeclarationElement(inputType, name, None) }
  }
}
