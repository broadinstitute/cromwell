package wdl.draft3.transforms.ast2wdlom

import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser
import wdl.draft3.parser.WdlParser.AstNode
import wdl.model.draft3.elements.OutputElement
import cats.syntax.apply._
import cats.syntax.either._
import common.transforms.CheckedAtoB

object CheckedAstToOutputElement {
  def convert(a: WdlParser.Ast): ErrorOr[OutputElement] = {
    val typeElementValidation: ErrorOr[String] = CheckedAtoB[AstNode, String].run(a.getAttribute("type")).toValidated
    val nameElementValidation: ErrorOr[String] = CheckedAtoB[AstNode, String].run(a.getAttribute("name")).toValidated
    val expressionElementValidation: ErrorOr[String] = CheckedAtoB[AstNode, String].run(a.getAttribute("expression")).toValidated

    (nameElementValidation, typeElementValidation, expressionElementValidation) mapN {
      (name, valueType, value) => OutputElement(valueType, name, value)
    }
  }
}
