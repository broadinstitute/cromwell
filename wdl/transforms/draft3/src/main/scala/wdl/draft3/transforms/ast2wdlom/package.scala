package wdl.draft3.transforms

import better.files.File
import cats.instances.either._
import common.transforms.CheckedAtoB
import common.validation.Checked._
import wdl.draft3.parser.WdlParser.{Ast, AstNode, Terminal}
import wdl.draft3.transforms.parsing._
import wdl.model.draft3.elements._

package object ast2wdlom {

  implicit val astNodeToAst: CheckedAtoB[AstNode, Ast] = CheckedAtoB.fromCheck {
    case ast: Ast => ast.validNelCheck
    case other => s"Cannot convert from AstNode type '${other.getClass.getSimpleName}' into Ast".invalidNelCheck
  }

  implicit val astToFileElement: CheckedAtoB[Ast, FileElement] = CheckedAtoB.fromErrorOr(CheckedAstToFileElement.convert)
  implicit val astNodeToImportElement: CheckedAtoB[AstNode, ImportElement] = astNodeToAst andThen CheckedAtoB.fromCheck(CheckedAstToImportElement.convert)
  implicit val astNodeToTaskDefinitionElement: CheckedAtoB[AstNode, TaskDefinitionElement] = astNodeToAst andThen CheckedAtoB.fromCheck(CheckedAstToTaskDefinitionElement.convert)
  implicit val astNodeToWorkflowDefinitionElement: CheckedAtoB[AstNode, WorkflowDefinitionElement] = astNodeToAst andThen CheckedAtoB.fromErrorOr(CheckedAstToWorkflowDefinitionElement.convert)
  implicit val checkedAstNodeToOutputsSectionElement: CheckedAtoB[AstNode, OutputsSectionElement] = astNodeToAst andThen CheckedAtoB.fromErrorOr(CheckedAstToOutputsSectionElement.convert)
  implicit val checkedAstNodeToOutputElement: CheckedAtoB[AstNode, OutputElement] = astNodeToAst andThen CheckedAtoB.fromErrorOr(CheckedAstToOutputElement.convert)

  implicit val fileToFileElement: CheckedAtoB[File, FileElement] = fileToAst andThen astToFileElement

  implicit val astNodeToString: CheckedAtoB[AstNode, String] = CheckedAtoB.fromCheck { (a: AstNode) => a match {
    case t: Terminal => t.getSourceString.validNelCheck
    case other: AstNode => s"Cannot convert ${other.getClass.getSimpleName} into String".invalidNelCheck
  }}
}
