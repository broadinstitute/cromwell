package wdl.draft3.transforms.ast2wdlom

import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser
import wdl.model.draft3.elements.{OutputElement, OutputsSectionElement}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import cats.syntax.either._
import cats.syntax.apply._
import wdl.draft3.parser.WdlParser.Ast


object CheckedAstToOutputsSectionElement {
  def convert(a: WdlParser.Ast): ErrorOr[OutputsSectionElement] = {
    val validatedWorkflowOutputDeclarations: ErrorOr[Vector[OutputElement]] =
      a.getAttributeAsVector[OutputElement]("outputs").toValidated

    validatedWorkflowOutputDeclarations map OutputsSectionElement

  }
}
