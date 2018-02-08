package wdl.draft3.transforms.ast2wdlom

import cats.syntax.apply._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements._
import wdl.draft3.parser.WdlParser.Ast
import wdl.draft3.transforms.ast2wdlom.FromAst._

object Draft3FileElementFromAst extends FromAst[FileElement] with FileElementMaker[Ast] {
  override def convert(ast: Ast): ErrorOr[FileElement] = if (ast.getName == "Draft3File") {

    val validatedImportElements: ErrorOr[Vector[ImportElement]] = ast.getAttributeAsVector[ImportElement]("imports")

    val validatedWorkflowAndTaskElements: ErrorOr[ (Vector[WorkflowDefinitionElement], Vector[TaskDefinitionElement]) ] =
      ast.getAttributeAsVectors[WorkflowDefinitionElement, TaskDefinitionElement]("body", "Workflow", "Task")

    (validatedImportElements, validatedWorkflowAndTaskElements) mapN { (imports, workflowAndTasks) =>
      FileElement(imports, workflowAndTasks._1, workflowAndTasks._2)
    }


  } else {
    s"Invalid AST type for FileElement: ${ast.getName}".invalidNel
  }
}
