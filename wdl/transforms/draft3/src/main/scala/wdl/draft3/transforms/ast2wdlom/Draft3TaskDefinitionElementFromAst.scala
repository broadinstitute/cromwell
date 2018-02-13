package wdl.draft3.transforms.ast2wdlom

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.Ast
import wdl.model.draft3.elements.{TaskDefinitionElement, TaskDefinitionElementMaker}

object Draft3TaskDefinitionElementFromAst extends FromAst[TaskDefinitionElement] with TaskDefinitionElementMaker[Ast] {
  override def convert(a: Ast): ErrorOr[TaskDefinitionElement] =
    "FromAst[TaskDefinitionElement](a: Ast) is not implemented".invalidNel
}
