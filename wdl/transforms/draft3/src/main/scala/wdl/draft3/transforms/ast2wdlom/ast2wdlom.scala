package wdl.draft3.transforms

import better.files.File
import cats.instances.either._
import common.transforms.CheckedAtoB
import common.validation.Checked._
import wdl.draft3.parser.WdlParser.{Ast, AstNode}
import wdl.model.draft3.elements.ExpressionElement.KvPair
import wdl.model.draft3.elements._
import wdl.transforms.base.ast2wdlom.{AstNodeToCommandPartElement, AstNodeToExpressionElement, AstNodeToKvPair, AstNodeToMetaKvPair, AstNodeToPlaceholderAttributeSet, AstNodeToStaticString, AstNodeToTypeElement, AstToCallElement, AstToCommandSectionElement, AstToDeclarationContent, AstToFileBodyElement, AstToFileElement, AstToIfElement, AstToImportElement, AstToInputDeclarationElement, AstToInputsSectionElement, AstToMetaSectionElement, AstToOutputsSectionElement, AstToParameterMetaSectionElement, AstToRuntimeAttributesSectionElement, AstToScatterElement, AstToStructElement, AstToTaskDefinitionElement, AstToTaskSectionElement, AstToWorkflowBodyElement, AstToWorkflowDefinitionElement, AstToWorkflowGraphNodeElementConverterMaker, GenericAst, GenericAstNode, astNodeToAst, astNodeToAstList}
import wdl.draft3.transforms.parsing.fileToAst
import wom.callable.MetaKvPair

package object ast2wdlom {

  val wrapAst: CheckedAtoB[Ast, GenericAst] = CheckedAtoB.fromCheck { a => Draft3GenericAst(a).validNelCheck }
  val wrapAstNode: CheckedAtoB[AstNode, GenericAstNode] = CheckedAtoB.fromCheck { a => Draft3GenericAstNode(a).validNelCheck }

  implicit val astNodeToStaticString: CheckedAtoB[GenericAstNode, StaticString] = AstNodeToStaticString.astNodeToStaticStringElement

  // meta sections
  implicit val astNodeToMetaKvPair: CheckedAtoB[GenericAstNode, MetaKvPair] = AstNodeToMetaKvPair.astNodeToMetaKvPair
  implicit val astNodeToMetaSectionElement: CheckedAtoB[GenericAstNode, MetaSectionElement] = astNodeToAst andThen AstToMetaSectionElement.astToMetaSectionElement
  implicit val astNodeToParameterMetaSectionElement: CheckedAtoB[GenericAstNode, ParameterMetaSectionElement] = astNodeToAst andThen AstToParameterMetaSectionElement.astToParameterMetaSectionElement

  implicit val astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement] = AstNodeToExpressionElement.astNodeToExpressionElement(customEngineFunctionMakers = Map.empty)
  implicit val astNodeToKvPair: CheckedAtoB[GenericAstNode, KvPair] = AstNodeToKvPair.astNodeToKvPair(astNodeToExpressionElement)

  implicit val astNodeToTypeElement: CheckedAtoB[GenericAstNode, TypeElement] = AstNodeToTypeElement.astNodeToTypeElement(Map.empty)
  implicit val astToStructElement: CheckedAtoB[GenericAst, StructElement] = AstToStructElement.astToStructElement
  implicit val astNodeToImportElement: CheckedAtoB[GenericAstNode, ImportElement] = astNodeToAst andThen AstToImportElement.astToImportElement

  implicit val astNodeToInputDeclarationElement: CheckedAtoB[GenericAstNode, InputDeclarationElement] = astNodeToAst andThen AstToInputDeclarationElement.astToInputDeclarationElement
  implicit val astNodeToInputsSectionElement: CheckedAtoB[GenericAstNode, InputsSectionElement] = astNodeToAst andThen AstToInputsSectionElement.astToInputsSectionElement

  implicit val astNodeToDeclarationContent: CheckedAtoB[GenericAstNode, DeclarationContent] = astNodeToAst andThen AstToDeclarationContent.astToDeclarationContent
  implicit val astNodeToOutputsSectionElement: CheckedAtoB[GenericAstNode, OutputsSectionElement] = astNodeToAst andThen AstToOutputsSectionElement.astToOutputSectionElement

  val astToWorkflowGraphNodeElementConverterMaker = new AstToWorkflowGraphNodeElementConverterMaker()
  implicit val astNodeToGraphElement: CheckedAtoB[GenericAstNode, WorkflowGraphElement] = astNodeToAst andThen astToWorkflowGraphNodeElementConverterMaker.converter
  implicit val astNodeToCallElement: CheckedAtoB[GenericAstNode, CallElement] = astNodeToAst andThen AstToCallElement.astToCallElement
  implicit val astNodeToScatterElement: CheckedAtoB[GenericAstNode, ScatterElement] = astNodeToAst andThen AstToScatterElement.astToScatterElement
  implicit val astNodeToIfElement: CheckedAtoB[GenericAstNode, IfElement] = astNodeToAst andThen AstToIfElement.astToIfElement
  astToWorkflowGraphNodeElementConverterMaker.astNodeToScatterElement = Some(astNodeToScatterElement)
  astToWorkflowGraphNodeElementConverterMaker.astNodeToIfElement = Some(astNodeToIfElement)
  astToWorkflowGraphNodeElementConverterMaker.astNodeToCallElement = Some(astNodeToCallElement)
  astToWorkflowGraphNodeElementConverterMaker.astNodeToDeclarationContent = Some(astNodeToDeclarationContent)

  implicit val astNodeToWorkflowBodyElement: CheckedAtoB[GenericAstNode, WorkflowBodyElement] = astNodeToAst andThen AstToWorkflowBodyElement.astToWorkflowBodyElement
  implicit val astToWorkflowDefinitionElement: CheckedAtoB[GenericAst, WorkflowDefinitionElement] = AstToWorkflowDefinitionElement.astToWorkflowDefinitionElement

  implicit val astNodeToPlaceholderAttributeSet: CheckedAtoB[GenericAstNode, PlaceholderAttributeSet] = astNodeToAstList andThen AstNodeToPlaceholderAttributeSet.attributeKvpConverter
  implicit val astNodeToCommandPartElement: CheckedAtoB[GenericAstNode, CommandPartElement] = AstNodeToCommandPartElement.astNodeToCommandPartElement
  implicit val astNodeToCommandSectionElement: CheckedAtoB[GenericAstNode, CommandSectionElement] = astNodeToAst andThen AstToCommandSectionElement.astToCommandSectionElement
  implicit val astNodeToRuntimeAttributesSectionElement: CheckedAtoB[GenericAstNode, RuntimeAttributesSectionElement] = astNodeToAst andThen AstToRuntimeAttributesSectionElement.astToRuntimeSectionElement
  implicit val astNodeToTaskSectionElement: CheckedAtoB[GenericAstNode, TaskSectionElement] = astNodeToAst andThen AstToTaskSectionElement.astToTaskSectionElement
  implicit val astToTaskDefinitionElement: CheckedAtoB[GenericAst, TaskDefinitionElement] = AstToTaskDefinitionElement.astToTaskDefinitionElement

  implicit val astToFileBodyElement: CheckedAtoB[GenericAstNode, FileBodyElement] = astNodeToAst andThen AstToFileBodyElement.astToFileBodyElement(astToWorkflowDefinitionElement, astToTaskDefinitionElement, astToStructElement)

  implicit val astToFileElement: CheckedAtoB[GenericAst, FileElement] = AstToFileElement.astToFileElement
  implicit val fileToFileElement: CheckedAtoB[File, FileElement] = fileToAst andThen wrapAst andThen astToFileElement

}
