package wdl.draft3.transforms

import better.files.File
import cats.instances.either._
import common.transforms.CheckedAtoB
import common.validation.Checked._
import wdl.draft3.parser.WdlParser.{Ast, AstNode, Terminal}
import wdl.draft3.transforms.parsing._
import wdl.model.draft3.elements.ExpressionElement.KvPair
import wdl.model.draft3.elements._

package object ast2wdlom {

  implicit val astNodeToAst: CheckedAtoB[AstNode, Ast] = CheckedAtoB.fromCheck {
    case ast: Ast => ast.validNelCheck
    case other => s"Cannot convert from AstNode type '${other.getClass.getSimpleName}' into Ast".invalidNelCheck
  }

  implicit val astToFileElement: CheckedAtoB[Ast, FileElement] = CheckedAtoB.fromErrorOr(AstToFileElement.convert)
  implicit val astToFileBodyElement: CheckedAtoB[AstNode, FileBodyElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToFileBodyElement.convert)
  implicit val astNodeToStructEntry: CheckedAtoB[AstNode, StructElement] = astNodeToAst andThen CheckedAtoB.fromErrorOr(AstToStructElement.convert)
  implicit val astNodeToImportElement: CheckedAtoB[AstNode, ImportElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToImportElement.convert)
  implicit val astNodeToTaskDefinitionElement: CheckedAtoB[AstNode, TaskDefinitionElement] = astNodeToAst andThen CheckedAtoB.fromErrorOr(AstToTaskDefinitionElement.convert)
  implicit val astNodeToWorkflowDefinitionElement: CheckedAtoB[AstNode, WorkflowDefinitionElement] = astNodeToAst andThen CheckedAtoB.fromErrorOr(AstToWorkflowDefinitionElement.convert)
  implicit val astNodeToInputsSectionElement: CheckedAtoB[AstNode, InputsSectionElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToInputsSectionElement.convert)
  implicit val astNodeToInputDeclarationElement: CheckedAtoB[AstNode, InputDeclarationElement] = astNodeToAst andThen CheckedAtoB.fromErrorOr(AstToInputDeclarationElement.convert)
  implicit val astNodeToTypeElement: CheckedAtoB[AstNode, TypeElement] = CheckedAtoB.fromErrorOr(AstNodeToTypeElement.convert)
  implicit val astNodeToExpressionElement: CheckedAtoB[AstNode, ExpressionElement] = CheckedAtoB.fromErrorOr(AstNodeToExpressionElement.convert)
  implicit val astNodeToKvPair: CheckedAtoB[AstNode, KvPair] = CheckedAtoB.fromErrorOr(AstNodeToKvPair.convert)
  implicit val astNodeToTaskSectionElement: CheckedAtoB[AstNode, TaskSectionElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToTaskSectionElement.convert)
  implicit val astNodeToWorkflowBodyElement: CheckedAtoB[AstNode, WorkflowBodyElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToWorkflowBodyElement.convert)
  implicit val astNodeToOutputsSectionElement: CheckedAtoB[AstNode, OutputsSectionElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToOutputsSectionElement.convert)
  implicit val astNodeToDeclarationContent: CheckedAtoB[AstNode, DeclarationContent] = astNodeToAst andThen CheckedAtoB.fromErrorOr(AstToDeclarationContent.convert)
  implicit val astNodeToRuntimeAttributesSectionElement: CheckedAtoB[AstNode, RuntimeAttributesSectionElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToRuntimeAttributesSectionElement.convert)
  implicit val astNodeToCommandSectionElement: CheckedAtoB[AstNode, CommandSectionElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToCommandSectionElement.convert)
  implicit val astNodeToCommandPartElement: CheckedAtoB[AstNode, CommandPartElement] = CheckedAtoB.fromErrorOr(AstNodeToCommandPartElement.convert)
  implicit val astNodeToCallElement: CheckedAtoB[AstNode, CallElement] = astNodeToAst andThen CheckedAtoB.fromErrorOr(AstToCallElement.convert)
  implicit val fileToFileElement: CheckedAtoB[File, FileElement] = fileToAst andThen astToFileElement

  implicit val astNodeToString: CheckedAtoB[AstNode, String] = CheckedAtoB.fromCheck { (a: AstNode) => a match {
    case t: Terminal => t.getSourceString.validNelCheck
    case other: AstNode => s"Cannot convert ${other.getClass.getSimpleName} into String".invalidNelCheck
  }}

  // meta sections
  implicit val astNodeToMetaValueElement: CheckedAtoB[AstNode, MetaValueElement] = CheckedAtoB.fromErrorOr(AstNodeToMetaValueElement.convert)
  implicit val astNodeToMetaKvPair: CheckedAtoB[AstNode, MetaKvPair] = CheckedAtoB.fromErrorOr(AstNodeToMetaKvPair.convert)
  implicit val astNodeToMetaSectionElement: CheckedAtoB[AstNode, MetaSectionElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToMetaSectionElement.convert)
  implicit val astNodeToParameterMetaSectionElement: CheckedAtoB[AstNode, ParameterMetaSectionElement] = astNodeToAst andThen CheckedAtoB.fromCheck(AstToParameterMetaSectionElement.convert)
}
