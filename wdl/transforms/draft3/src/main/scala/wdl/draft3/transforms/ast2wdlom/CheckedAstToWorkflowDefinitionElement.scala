package wdl.draft3.transforms.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.{Ast, AstNode}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.{WorkflowDefinitionElement, OutputsSectionElement}

object CheckedAstToWorkflowDefinitionElement {

  def convert(a: Ast): ErrorOr[WorkflowDefinitionElement] = {

    val outputsSectionElementValidation: ErrorOr[Vector[OutputsSectionElement]] = a.getAttributeAsVector[OutputsSectionElement]("body").toValidated

    val nameElementValidation: ErrorOr[String] = CheckedAtoB[AstNode, String].run(a.getAttribute("name")).toValidated

    (nameElementValidation, outputsSectionElementValidation) mapN { (name, outputs) => WorkflowDefinitionElement(name, outputs) }
  }
}
