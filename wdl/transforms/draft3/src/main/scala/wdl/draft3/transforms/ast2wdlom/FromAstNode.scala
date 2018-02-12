package wdl.draft3.transforms.ast2wdlom

import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.AstNode


trait FromAstNode[A] extends FromAtoB[AstNode, A]

object FromAstNode {
  def apply[A](ast: AstNode)(implicit value: FromAtoB[AstNode, A]): ErrorOr[A] = { value.convert(ast) }
}

