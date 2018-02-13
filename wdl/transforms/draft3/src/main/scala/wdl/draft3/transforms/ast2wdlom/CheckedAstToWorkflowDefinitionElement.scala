package wdl.draft3.transforms.ast2wdlom

import cats.data.NonEmptyList
import cats.syntax.apply._
import cats.syntax.either._
import common.Checked
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.{Ast, AstNode}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.WorkflowDefinitionElement

object CheckedAstToWorkflowDefinitionElement {

  type CheckedAstToWorkflowDefinitionElement = CheckedAstTo[WorkflowDefinitionElement]
  def instance: CheckedAstToWorkflowDefinitionElement = CheckedAtoB(convert _)

  def convert(a: Ast): Checked[WorkflowDefinitionElement] = {

    val noUnexpectedElementsValidation: ErrorOr[Unit] = (a.getAttributeAsAstNodeVector("body") flatMap {
      case e if e.isEmpty => Right(())
      case nel if nel.nonEmpty => Left(NonEmptyList.fromListUnsafe(nel.toList).map(ast => s"Unexpected item in Workflow AST: ${ast.getClass.getSimpleName}"))
    }).toValidated

    val nameElementValidation: ErrorOr[String] = CheckedAtoB[AstNode, String].run(a.getAttribute("name")).toValidated

    val resultErrorOr = (nameElementValidation, noUnexpectedElementsValidation) mapN { (name, _) => WorkflowDefinitionElement(name) }
    resultErrorOr.toEither
  }
}
