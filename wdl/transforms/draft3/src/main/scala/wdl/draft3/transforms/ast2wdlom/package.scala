package wdl.draft3.transforms

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.{AstNode, Terminal}
import wdl.draft3.transforms.parsing.FileParser

package object ast2wdlom {

  implicit val astFromAstNode = AstFromAstNode
  implicit val draft3FileElementFromAstNode = FromAtoB.viaX(AstFromAstNode, Draft3FileElementFromAst)
  implicit val draft3ImportElementFromAstNode = FromAtoB.viaX(AstFromAstNode, Draft3ImportElementFromAst)
  implicit val draft3TaskDefinitionElementFromAstNode = FromAtoB.viaX(AstFromAstNode, Draft3TaskDefinitionElementFromAst)
  implicit val draft3WorkflowDefinitionElementFromAstNode = FromAtoB.viaX(AstFromAstNode, Draft3WorkflowDefinitionElementFromAst)

  implicit val draft3FileElementFromFile = FromAtoB.viaX(FileParser, Draft3FileElementFromAst)

  implicit val StringFromAstNode: FromAtoB[AstNode, String] = new FromAtoB[AstNode, String] {
    override def convert(a: AstNode): ErrorOr[String] = a match {
      case t: Terminal => t.getSourceString.valid
      case other: AstNode => s"Cannot convert ${other.getClass.getSimpleName} into String".invalidNel
    }
  }
}
