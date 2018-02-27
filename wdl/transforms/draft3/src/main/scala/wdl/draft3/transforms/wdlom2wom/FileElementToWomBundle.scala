package wdl.draft3.transforms.wdlom2wom

import cats.instances.vector._
import cats.syntax.either._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.{ErrorOr, _}
import wdl.model.draft3.elements.{FileElement, ImportElement, TaskDefinitionElement, WorkflowDefinitionElement}
import wom.callable.{TaskDefinition, WorkflowDefinition}
import wom.executable.WomBundle
import wom.transforms.WomBundleMaker
import wom.transforms.WomBundleMaker.ops._

import scala.concurrent.Future

object FileElementToWomBundle {

  implicit val fileElementToWomBundle: WomBundleMaker[FileElement] = new WomBundleMaker[FileElement] {

    override def toWomBundle(a: FileElement, importResolvers: List[String => Future[Checked[WomBundle]]]): Checked[WomBundle] = {

      val importsValidation: ErrorOr[Vector[ImportElement]] = if (a.imports.isEmpty) Vector.empty.valid else "FileElement to WOM conversion of imports not yet implemented.".invalidNel
      val tasksValidation: ErrorOr[Vector[TaskDefinitionElement]] = if (a.imports.isEmpty) Vector.empty.valid else "FileElement to WOM conversion of tasks not yet implemented.".invalidNel


      def importAndTasksToWorkflow(imports: Vector[ImportElement], tasks: Vector[TaskDefinitionElement]): ErrorOr[WomBundle] = {
        implicit val workflowConverter: CheckedAtoB[WorkflowDefinitionElement, WorkflowDefinition] = workflowDefinitionElementToWomWorkflowDefinition

        val tasks: Set[TaskDefinition] = Set.empty

        val workflowsValidation: ErrorOr[Vector[WorkflowDefinition]] = {
          a.workflows.toVector.traverse[ErrorOr, WorkflowDefinition](workflowConverter.run(_).toValidated)
        }

        workflowsValidation map { workflows =>
          WomBundle(tasks ++ workflows, Map.empty)
        }
      }

      ((importsValidation, tasksValidation) flatMapN { importAndTasksToWorkflow }).toEither
    }
  }

  def convert(a: FileElementAndImportResolvers): Checked[WomBundle] = a.fileElement.toWomBundle(a.importResolvers)
}

final case class FileElementAndImportResolvers(fileElement: FileElement, importResolvers: List[String => Future[Checked[WomBundle]]])