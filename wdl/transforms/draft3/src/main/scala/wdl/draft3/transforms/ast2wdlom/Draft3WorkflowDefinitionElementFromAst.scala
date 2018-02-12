package wdl.draft3.transforms.ast2wdlom

import cats.data.NonEmptyList
import cats.syntax.apply._
import cats.syntax.either._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.{Ast, AstNode}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.{WorkflowDefinitionElement, WorkflowDefinitionElementMaker}

object Draft3WorkflowDefinitionElementFromAst extends FromAst[WorkflowDefinitionElement] with WorkflowDefinitionElementMaker[Ast] {
  override def convert(a: Ast): ErrorOr[WorkflowDefinitionElement] = {

    val noUnexpectedElementsValidation: ErrorOr[Unit] = (a.getAttributeAsAstNodeVector("body") flatMap {
      case e if e.isEmpty => Right(())
      case nel if nel.nonEmpty => Left(NonEmptyList.fromListUnsafe(nel.toList).map(ast => s"Unexpected item in Workflow AST: ${ast.getClass.getSimpleName}"))
    }).toValidated

    val nameElementValidation: ErrorOr[String] = FromAtoB[AstNode, String].convert(a.getAttribute("name"))

    (nameElementValidation, noUnexpectedElementsValidation) mapN { (name, _) => WorkflowDefinitionElement(name) }

  }
}
