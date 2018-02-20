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

  implicit val astToFileElement: CheckedAtoB[Ast, FileElement] = CheckedAtoB.fromErrorOr(AstToFileElement.convert)
  implicit val astToFileBodyElement: CheckedAtoB[AstNode, FileBodyElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToFileBodyElement.convert)
  implicit val astNodeToImportElement: CheckedAtoB[AstNode, ImportElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToImportElement.convert)
  implicit val astNodeToTaskDefinitionElement: CheckedAtoB[AstNode, TaskDefinitionElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToTaskDefinitionElement.convert)
  implicit val astNodeToWorkflowDefinitionElement: CheckedAtoB[AstNode, WorkflowDefinitionElement] = astNodeToAst andThen CheckedAtoB.fromErrorOr(AstToWorkflowDefinitionElement.convert)
  implicit val astNodeToInputsSectionElement: CheckedAtoB[AstNode, InputsSectionElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToInputsSectionElement.convert)
  implicit val astNodeToInputDeclarationElement: CheckedAtoB[AstNode, InputDeclarationElement] = astNodeToAst andThen CheckedAtoB.fromErrorOr(AstToInputDeclarationElement.convert)
  implicit val astNodeToTypeElement: CheckedAtoB[AstNode, TypeElement] = CheckedAtoB.fromErrorOr(AstNodeToTypeElement.convert)
  implicit val astNodeToExpressionElement: CheckedAtoB[AstNode, ExpressionElement] = CheckedAtoB.fromErrorOr(AstNodeToExpressionElement.convert)
  implicit val astNodeToTaskBodyElement: CheckedAtoB[AstNode, TaskBodyElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToTaskBodyElement.convert)
  implicit val astNodeToWorkflowBodyElement: CheckedAtoB[AstNode, WorkflowBodyElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToWorkflowBodyElement.convert)
  implicit val astNodeToOutputsSectionElement: CheckedAtoB[AstNode, OutputsSectionElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToOutputsSectionElement.convert)
  implicit val astNodeToDeclarationContent: CheckedAtoB[AstNode, DeclarationContent] = astNodeToAst andThen CheckedAtoB.fromErrorOr(AstToDeclarationContent.convert)

  implicit val fileToFileElement: CheckedAtoB[File, FileElement] = fileToAst andThen astToFileElement

  implicit val astNodeToString: CheckedAtoB[AstNode, String] = CheckedAtoB.fromCheck { (a: AstNode) => a match {
    case t: Terminal => t.getSourceString.validNelCheck
    case other: AstNode => s"Cannot convert ${other.getClass.getSimpleName} into String".invalidNelCheck
  }}
}
