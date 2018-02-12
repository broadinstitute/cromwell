package wdl.draft3.transforms.wdlom2wom

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import wdl.draft3.transforms.ast2wdlom.FromAtoB
import wdl.model.draft3.elements.{FileElement, ImportElement, TaskDefinitionElement, WorkflowDefinitionElement}
import wom.callable.WorkflowDefinition
import wom.executable.Executable

case class FileElementToWomExecutable(inputs: Option[String]) extends FromAtoB[FileElement, Executable] {
  override def convert(a: FileElement): ErrorOr[Executable] = {

    val importsValidation: ErrorOr[Vector[ImportElement]] = if (a.imports.isEmpty) Vector.empty.valid else "FileElement to WOM conversion of imports not yet implemented.".invalidNel
    val tasksValidation: ErrorOr[Vector[TaskDefinitionElement]] = if (a.imports.isEmpty) Vector.empty.valid else "FileElement to WOM conversion of tasks not yet implemented.".invalidNel


    (importsValidation, tasksValidation) flatMapN { (_, _) =>
      implicit val workflowConverter: FromAtoB[WorkflowDefinitionElement, WorkflowDefinition] = WorkflowDefinitionElementToWomWorkflowDefinition(inputs)
      val workflowsValidation: ErrorOr[Vector[WorkflowDefinition]] = FromAtoB.forVectors[WorkflowDefinitionElement, WorkflowDefinition].convert(a.workflows.toVector)

      workflowsValidation flatMap {
        case one if one.size == 1 => Executable(one.head, Map.empty).valid
        case notOne => s"Cannot turn a ${notOne.size}-workflow WDL FileElement into a WomExecutable (must be exactly 1)".invalidNel
      }
    }
  }
}
