package wdl.transforms.base.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr._
import wdl.model.draft3.elements._

object AstToFileElement {

  def astToFileElement(implicit astNodeToImportElement: CheckedAtoB[GenericAstNode, ImportElement],
                       astNodeToFileBodyElement: CheckedAtoB[GenericAstNode, FileBodyElement]
                      ): CheckedAtoB[GenericAst, FileElement] = CheckedAtoB.fromErrorOr { ast =>

    val validatedImportElements: ErrorOr[Vector[ImportElement]] = ast.getAttributeAsVector[ImportElement]("imports").toValidated
    val validatedFileBodyElements: ErrorOr[Vector[FileBodyElement]] = ast.getAttributeAsVector[FileBodyElement]("body").toValidated

    (validatedImportElements, validatedFileBodyElements) mapN { (importElements, fileBodyElements) =>
      val workflowElements: Vector[WorkflowDefinitionElement] = fileBodyElements.collect { case e: WorkflowDefinitionElement => e }
      val taskElements: Vector[TaskDefinitionElement] = fileBodyElements.collect { case e: TaskDefinitionElement => e }
      val structElements: Vector[StructElement] = fileBodyElements.collect { case e: StructElement => e }
      FileElement(importElements, structElements, workflowElements, taskElements)
    }
  }
}
