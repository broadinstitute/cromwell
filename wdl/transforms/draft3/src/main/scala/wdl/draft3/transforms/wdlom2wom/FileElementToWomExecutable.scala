package wdl.draft3.transforms.wdlom2wom

import cats.syntax.validated._
import cats.syntax.either._
import cats.syntax.traverse._
import cats.instances.vector._
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import common.transforms.CheckedAtoB
import wdl.model.draft3.elements.{FileElement, ImportElement, TaskDefinitionElement, WorkflowDefinitionElement}
import wom.callable.WorkflowDefinition
import wom.executable.Executable
import wom.transforms.WomExecutableMaker.ExecutableMakerInputs

object FileElementToWomExecutable {

  def convert(a: ExecutableMakerInputs[FileElement]): ErrorOr[Executable] = {

    val importsValidation: ErrorOr[Vector[ImportElement]] = if (a.from.imports.isEmpty) Vector.empty.valid else "FileElement to WOM conversion of imports not yet implemented.".invalidNel
    val tasksValidation: ErrorOr[Vector[TaskDefinitionElement]] = if (a.from.imports.isEmpty) Vector.empty.valid else "FileElement to WOM conversion of tasks not yet implemented.".invalidNel


    def importAndTasksToWorkflow(imports: Vector[ImportElement], tasks: Vector[TaskDefinitionElement]): ErrorOr[Executable] = {
      implicit val workflowConverter: CheckedAtoB[WorkflowDefinitionElement, WorkflowDefinition] = workflowDefinitionElementToWomWorkflowDefinition

      val workflowsValidation: ErrorOr[Vector[WorkflowDefinition]] = {
        a.from.workflows.toVector.traverse[ErrorOr, WorkflowDefinition](workflowConverter.run(_).toValidated)
      }

      workflowsValidation flatMap {
        case one if one.size == 1 => Executable(one.head, Map.empty).validNel
        case notOne => s"Cannot turn a ${notOne.size}-workflow WDL FileElement into a WomExecutable (must be exactly 1)".invalidNel
      }
    }


    (importsValidation, tasksValidation) flatMapN importAndTasksToWorkflow

  }
}
