package wdl.draft3.transforms.wdlom2wom

import cats.syntax.validated._
import cats.syntax.either._
import cats.syntax.traverse._
import cats.instances.vector._
import common.Checked
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import wdl.draft3.transforms.ast2wdlom.CheckedAtoB
import wdl.model.draft3.elements.{FileElement, ImportElement, TaskDefinitionElement, WorkflowDefinitionElement}
import wom.callable.WorkflowDefinition
import wom.executable.Executable

object FileElementToWomExecutable {

  type FileElementToWomExecutable = CheckedAtoB[FileElementAndInputsFile, Executable]
  def instance: FileElementToWomExecutable = CheckedAtoB(convert _)

  final case class FileElementAndInputsFile(fileElement: FileElement, inputs: Option[String])

  def convert(a: FileElementAndInputsFile): Checked[Executable] = {

    val importsValidation: ErrorOr[Vector[ImportElement]] = if (a.fileElement.imports.isEmpty) Vector.empty.valid else "FileElement to WOM conversion of imports not yet implemented.".invalidNel
    val tasksValidation: ErrorOr[Vector[TaskDefinitionElement]] = if (a.fileElement.imports.isEmpty) Vector.empty.valid else "FileElement to WOM conversion of tasks not yet implemented.".invalidNel


    def importAndTasksToWorkflow(imports: Vector[ImportElement], tasks: Vector[TaskDefinitionElement]): ErrorOr[Executable] = {
      implicit val workflowConverter: CheckedAtoB[WorkflowDefinitionElement, WorkflowDefinition] = WorkflowDefinitionElementToWomWorkflowDefinition.instance

      val workflowsValidation: ErrorOr[Vector[WorkflowDefinition]] = {
        a.fileElement.workflows.toVector.traverse[ErrorOr, WorkflowDefinition](workflowConverter.run(_).toValidated)
      }

      workflowsValidation flatMap {
        case one if one.size == 1 => Executable(one.head, Map.empty).validNel
        case notOne => s"Cannot turn a ${notOne.size}-workflow WDL FileElement into a WomExecutable (must be exactly 1)".invalidNel
      }
    }


    ((importsValidation, tasksValidation) flatMapN importAndTasksToWorkflow).toEither

  }
}
