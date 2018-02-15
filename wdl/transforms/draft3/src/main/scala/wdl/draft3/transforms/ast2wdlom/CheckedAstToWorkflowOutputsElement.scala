package wdl.draft3.transforms.ast2wdlom

import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser
import wdl.model.draft3.elements.{WorkflowOutputDeclarationElement, WorkflowOutputsElement}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import cats.syntax.either._
import cats.syntax.apply._
import wdl.draft3.parser.WdlParser.Ast


object CheckedAstToWorkflowOutputsElement {
  def convert(a: WdlParser.Ast): ErrorOr[WorkflowOutputsElement] = {
    val validatedWorkflowOutputDeclarations: ErrorOr[Vector[WorkflowOutputDeclarationElement]] =
      a.getAttributeAsVector[WorkflowOutputDeclarationElement]("outputs").toValidated

    validatedWorkflowOutputDeclarations map WorkflowOutputsElement

  }
}
