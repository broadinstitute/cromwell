package wdl.draft3.transforms.ast2wdlom

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.Ast
import wdl.model.draft3.elements.{WorkflowDefinitionElement, WorkflowDefinitionElementMaker}

object Draft3WorkflowDefinitionElementFromAst extends FromAst[WorkflowDefinitionElement] with WorkflowDefinitionElementMaker[Ast] {
  override def convert(a: Ast): ErrorOr[WorkflowDefinitionElement] =
    "FromAst[WorkflowDefinitionElement](a: Ast) is not implemented".invalidNel
}
