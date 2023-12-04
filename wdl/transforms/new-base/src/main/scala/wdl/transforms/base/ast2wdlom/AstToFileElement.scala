package wdl.transforms.base.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import common.collections.EnhancedCollections._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr._
import wdl.model.draft3.elements._

object AstToFileElement {

  def astToFileElement(implicit
    astNodeToImportElement: CheckedAtoB[GenericAstNode, ImportElement],
    astNodeToFileBodyElement: CheckedAtoB[GenericAstNode, FileBodyElement]
  ): CheckedAtoB[GenericAst, FileElement] = CheckedAtoB.fromErrorOr { ast =>
    val validatedImportElements: ErrorOr[Vector[ImportElement]] =
      ast.getAttributeAsVector[ImportElement]("imports").toValidated
    val validatedFileBodyElements: ErrorOr[Vector[FileBodyElement]] =
      ast.getAttributeAsVector[FileBodyElement]("body").toValidated

    (validatedImportElements, validatedFileBodyElements) mapN { (importElements, fileBodyElements) =>
      val workflowElements: Vector[WorkflowDefinitionElement] = fileBodyElements.filterByType[WorkflowDefinitionElement]
      val taskElements: Vector[TaskDefinitionElement] = fileBodyElements.filterByType[TaskDefinitionElement]
      val structElements: Vector[StructElement] = fileBodyElements.filterByType[StructElement]
      FileElement(importElements, structElements, workflowElements, taskElements)
    }
  }
}
