package wdl.draft3.transforms.ast2wdlom

import common.Checked
import common.validation.Checked._
import wdl.draft3.parser.WdlParser.{Ast, AstNode}

object CheckedAstNodeToAst {

  type CheckedAstNodeToAst = CheckedAstNodeTo[Ast]
  def instance: CheckedAstNodeToAst = CheckedAtoB(convert _)

  def convert(a: AstNode): Checked[Ast] = a match {
    case ast: Ast => ast.validNelCheck
    case other => s"Cannot convert from AstNode type '${other.getClass.getSimpleName}' into Ast".invalidNelCheck
  }
}
