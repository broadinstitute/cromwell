package wdl.transforms.v1_2

import cats.instances.either._
import better.files.File
import common.transforms.CheckedAtoB
import common.validation.Checked._
import wdl.v1_2.parser.WdlParser.{Ast, AstNode}
import wdl.model.draft3.elements.ExpressionElement.KvPair
import wdl.model.draft3.elements._
import wdl.transforms.base.ast2wdlom._
import wdl.transforms.v1_2.parsing.fileToAst
import wom.callable.MetaKvPair
import wom.types.WomUnlistedDirectoryType

package object ast2wdlom {

  val wrapAst: CheckedAtoB[Ast, GenericAst] = CheckedAtoB.fromCheck(a => V1_2GenericAst(a).validNelCheck)
  val wrapAstNode: CheckedAtoB[AstNode, GenericAstNode] = CheckedAtoB.fromCheck { a =>
    V1_2GenericAstNode(a).validNelCheck
  }

  implicit val astNodeToStaticString: CheckedAtoB[GenericAstNode, StaticString] =
    AstNodeToStaticString.astNodeToStaticStringElement()

  // meta sections
  implicit val astNodeToMetaKvPair: CheckedAtoB[GenericAstNode, MetaKvPair] = AstNodeToMetaKvPair.astNodeToMetaKvPair
  implicit val astNodeToMetaSectionElement: CheckedAtoB[GenericAstNode, MetaSectionElement] =
    astNodeToAst andThen AstToMetaSectionElement.astToMetaSectionElement
  implicit val astNodeToParameterMetaSectionElement: CheckedAtoB[GenericAstNode, ParameterMetaSectionElement] =
    astNodeToAst andThen AstToParameterMetaSectionElement.astToParameterMetaSectionElement

  implicit val astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement] =
    AstNodeToExpressionElement.astNodeToExpressionElement(customEngineFunctionMakers =
      AstToNewExpressionElements.newV1_2EngineFunctionMakers
    )
  implicit val astNodeToKvPair: CheckedAtoB[GenericAstNode, KvPair] =
    V1_2AstNodeToKvPair.astNodeToKvPair(astNodeToExpressionElement)

  // WDL 1.2 adds the `Directory` type
  implicit val astNodeToTypeElement: CheckedAtoB[GenericAstNode, TypeElement] =
    AstNodeToTypeElement.astNodeToTypeElement(Map("Directory" -> WomUnlistedDirectoryType))
  // WDL 1.2 structs may have `meta` and `parameter_meta` sections
  implicit val astToStructElement: CheckedAtoB[GenericAst, StructElement] = AstToV1_2Elements.astToStructElement
  implicit val astNodeToImportElement: CheckedAtoB[GenericAstNode, ImportElement] =
    astNodeToAst andThen AstToImportElement.astToImportElement

  implicit val astNodeToInputDeclarationElement: CheckedAtoB[GenericAstNode, InputDeclarationElement] =
    astNodeToAst andThen AstToInputDeclarationElement.astToInputDeclarationElement
  implicit val astNodeToInputsSectionElement: CheckedAtoB[GenericAstNode, InputsSectionElement] =
    astNodeToAst andThen AstToInputsSectionElement.astToInputsSectionElement

  implicit val astNodeToDeclarationContent: CheckedAtoB[GenericAstNode, DeclarationContent] =
    astNodeToAst andThen AstToDeclarationContent.astToDeclarationContent
  implicit val astNodeToOutputsSectionElement: CheckedAtoB[GenericAstNode, OutputsSectionElement] =
    astNodeToAst andThen AstToOutputsSectionElement.astToOutputSectionElement

  val astToWorkflowGraphNodeElementConverterMaker = new AstToWorkflowGraphNodeElementConverterMaker()
  implicit val astNodeToGraphElement: CheckedAtoB[GenericAstNode, WorkflowGraphElement] =
    astNodeToAst andThen astToWorkflowGraphNodeElementConverterMaker.converter
  implicit val astNodeToCallElement: CheckedAtoB[GenericAstNode, CallElement] =
    astNodeToAst andThen AstToCallElement.astToCallElement
  implicit val astNodeToScatterElement: CheckedAtoB[GenericAstNode, ScatterElement] =
    astNodeToAst andThen AstToScatterElement.astToScatterElement
  implicit val astNodeToIfElement: CheckedAtoB[GenericAstNode, IfElement] =
    astNodeToAst andThen AstToIfElement.astToIfElement
  astToWorkflowGraphNodeElementConverterMaker.astNodeToScatterElement = Some(astNodeToScatterElement)
  astToWorkflowGraphNodeElementConverterMaker.astNodeToIfElement = Some(astNodeToIfElement)
  astToWorkflowGraphNodeElementConverterMaker.astNodeToCallElement = Some(astNodeToCallElement)
  astToWorkflowGraphNodeElementConverterMaker.astNodeToDeclarationContent = Some(astNodeToDeclarationContent)

  // WDL 1.2 adds the workflow `hints` section
  implicit val astNodeToWorkflowBodyElement: CheckedAtoB[GenericAstNode, WorkflowBodyElement] =
    astNodeToAst andThen AstToV1_2Elements.astToWorkflowBodyElement(AstToWorkflowBodyElement.astToWorkflowBodyElement)
  implicit val astToWorkflowDefinitionElement: CheckedAtoB[GenericAst, WorkflowDefinitionElement] =
    AstToV1_2Elements.astToWorkflowDefinitionElement(AstToWorkflowDefinitionElement.astToWorkflowDefinitionElement)

  implicit val astNodeToPlaceholderAttributeSet: CheckedAtoB[GenericAstNode, PlaceholderAttributeSet] =
    astNodeToAstList andThen AstNodeToPlaceholderAttributeSet.attributeKvpConverter
  implicit val astNodeToCommandPartElement: CheckedAtoB[GenericAstNode, CommandPartElement] =
    AstNodeToCommandPartElement.astNodeToCommandPartElement
  implicit val astNodeToCommandSectionElement: CheckedAtoB[GenericAstNode, CommandSectionElement] =
    astNodeToAst andThen AstToCommandSectionElement.astToCommandSectionElement
  implicit val astNodeToRuntimeAttributesSectionElement: CheckedAtoB[GenericAstNode, RuntimeAttributesSectionElement] =
    astNodeToAst andThen AstToRuntimeAttributesSectionElement.astToRuntimeSectionElement
  // WDL 1.2 adds the task `requirements` and `hints` sections and the `env` declaration modifier
  implicit val astNodeToTaskSectionElement: CheckedAtoB[GenericAstNode, TaskSectionElement] =
    astNodeToAst andThen AstToV1_2Elements.astToTaskSectionElement(AstToTaskSectionElement.astToTaskSectionElement)
  implicit val astToTaskDefinitionElement: CheckedAtoB[GenericAst, TaskDefinitionElement] =
    AstToV1_2Elements.astToTaskDefinitionElement

  implicit val astToFileBodyElement: CheckedAtoB[GenericAstNode, FileBodyElement] =
    astNodeToAst andThen AstToFileBodyElement.astToFileBodyElement(astToWorkflowDefinitionElement,
                                                                   astToTaskDefinitionElement,
                                                                   astToStructElement
    )

  implicit val astToFileElement: CheckedAtoB[GenericAst, FileElement] = AstToFileElement.astToFileElement
  implicit val fileToFileElement: CheckedAtoB[File, FileElement] = fileToAst andThen wrapAst andThen astToFileElement

}
