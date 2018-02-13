package wdl.draft3.transforms.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import common.Checked
import common.validation.ErrorOr.ErrorOr
import common.validation.Checked._
import wdl.model.draft3.elements._
import wdl.draft3.parser.WdlParser.Ast
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._

object CheckedAstToFileElement {

  type CheckedAstToFileElement = CheckedAstTo[FileElement]
  def instance: CheckedAstToFileElement = CheckedAtoB(convert _)

  def convert(ast: Ast): Checked[FileElement] = if (ast.getName == "Draft3File") {

    val validatedImportElements: ErrorOr[Vector[ImportElement]] = ast.getAttributeAsVector[ImportElement]("imports").toValidated

    val validatedWorkflowAndTaskElements: ErrorOr[ (Vector[WorkflowDefinitionElement], Vector[TaskDefinitionElement]) ] =
      ast.getAttributeAsVectors[WorkflowDefinitionElement, TaskDefinitionElement]("body", "Workflow", "Task").toValidated

    val result = (validatedImportElements, validatedWorkflowAndTaskElements) mapN { (imports, workflowAndTasks) =>
      FileElement(imports, workflowAndTasks._1, workflowAndTasks._2)
    }

    result.toEither

  } else {
    s"Invalid AST type for FileElement: ${ast.getName}".invalidNelCheck
  }
}
