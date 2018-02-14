package wdl.draft3.transforms

import better.files.File
import cats.instances.either._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.Checked._
import wdl.draft3.parser.WdlParser.{Ast, AstNode, Terminal}
import wdl.draft3.transforms.parsing._
import wdl.model.draft3.elements.FileElement

package object ast2wdlom {

  implicit val checkedAstNodeToAst = CheckedAtoB.fromCheck(run = (a: AstNode) => a match {
    case ast: Ast => ast.validNelCheck
    case other => s"Cannot convert from AstNode type '${other.getClass.getSimpleName}' into Ast".invalidNelCheck
  })

  implicit val checkedAstToFileElement = CheckedAtoB.fromErrorOr(CheckedAstToFileElement.convert)
  implicit val checkedAstNodeToImportElement = checkedAstNodeToAst andThen CheckedAtoB.fromCheck(CheckedAstToImportElement.convert)
  implicit val checkedAstNodeToTaskDefinitionElement = checkedAstNodeToAst andThen CheckedAtoB.fromCheck(CheckedAstToTaskDefinitionElement.convert)
  implicit val checkedAstNodeToWorkflowDefinitionElement = checkedAstNodeToAst andThen CheckedAtoB.fromErrorOr(CheckedAstToWorkflowDefinitionElement.convert)

  implicit val checkedFileElementFromFile: CheckedAtoB[File, FileElement] = checkedFileToAst andThen checkedAstToFileElement

  implicit val checkedAstNodeToString = CheckedAtoB.fromCheck { (a: AstNode) => a match {
    case t: Terminal => t.getSourceString.validNelCheck
    case other: AstNode => s"Cannot convert ${other.getClass.getSimpleName} into String".invalidNelCheck
  }}

}
