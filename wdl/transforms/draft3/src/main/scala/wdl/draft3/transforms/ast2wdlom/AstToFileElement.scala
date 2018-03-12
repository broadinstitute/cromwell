package wdl.draft3.transforms.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.collections.EnhancedCollections._
import wdl.model.draft3.elements._
import wdl.draft3.parser.WdlParser.Ast
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._

object AstToFileElement {

  def convert(ast: Ast): ErrorOr[FileElement] = if (ast.getName == "Draft3File") {

    val validatedImportElements: ErrorOr[Vector[ImportElement]] = ast.getAttributeAsVector[ImportElement]("imports").toValidated
    val validatedFileBodyElements: ErrorOr[Vector[FileBodyElement]] = ast.getAttributeAsVector[FileBodyElement]("body").toValidated

    (validatedImportElements, validatedFileBodyElements) mapN { (importElements, fileBodyElements) =>
      val workflowElements: Vector[WorkflowDefinitionElement] = fileBodyElements.filterByType[WorkflowDefinitionElement]
      val taskElements: Vector[TaskDefinitionElement] = fileBodyElements.filterByType[TaskDefinitionElement]
      val structElements: Vector[StructElement] = fileBodyElements.filterByType[StructElement]
      FileElement(importElements, structElements, workflowElements, taskElements)
    }

  } else {
    s"Invalid AST type for FileElement: ${ast.getName}".invalidNel
  }
}
