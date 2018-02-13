package wdl.draft3.transforms.ast2wdlom

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.{Ast, AstNode}

object AstFromAstNode extends FromAstNode[Ast] {
  override def convert(a: AstNode): ErrorOr[Ast] = a match {
    case ast: Ast => ast.valid
    case other => s"Cannot convert from AstNode type '${other.getClass.getSimpleName}' into Ast".invalidNel
  }
}
