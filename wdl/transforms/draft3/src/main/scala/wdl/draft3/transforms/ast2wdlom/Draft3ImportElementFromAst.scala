package wdl.draft3.transforms.ast2wdlom

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.Ast
import wdl.model.draft3.elements.{ImportElement, ImportElementMaker}

object Draft3ImportElementFromAst extends FromAst[ImportElement] with ImportElementMaker[Ast] {
  override def convert(a: Ast): ErrorOr[ImportElement] =
    "FromAst[ImportElement](a: Ast) is not implemented".invalidNel
}
