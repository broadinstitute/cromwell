package wdl.draft3.transforms.ast2wdlom

import cats.data.NonEmptyList
import cats.syntax.apply._
import cats.syntax.either._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.{Ast, AstNode}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.{WorkflowDefinitionElement, WorkflowOutputsElement}

object CheckedAstToWorkflowDefinitionElement {

  def convert(a: Ast): ErrorOr[WorkflowDefinitionElement] = {

    val workflowOutputsElementValidation: ErrorOr[Vector[WorkflowOutputsElement]] = a.getAttributeAsVector[WorkflowOutputsElement]("body").toValidated

    val nameElementValidation: ErrorOr[String] = CheckedAtoB[AstNode, String].run(a.getAttribute("name")).toValidated

    (nameElementValidation, workflowOutputsElementValidation) mapN { (name, workflowOutputs) => WorkflowDefinitionElement(name, workflowOutputs) }
  }
}
