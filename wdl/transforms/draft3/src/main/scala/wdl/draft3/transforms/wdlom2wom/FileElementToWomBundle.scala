package wdl.draft3.transforms.wdlom2wom

import cats.instances.vector._
import cats.syntax.either._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import wdl.draft3.transforms.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition.WorkflowDefinitionConvertInputs
import wdl.model.draft3.elements._
import wom.callable.{TaskDefinition, WorkflowDefinition}
import wom.executable.WomBundle
import wom.transforms.WomBundleMaker
import wom.transforms.WomBundleMaker.ops._
import wdl.draft3.transforms.wdlom2wom.StructEvaluation.StructEvaluationInputs
import wom.types.WomType

import scala.concurrent.Future

object FileElementToWomBundle {

  implicit val fileElementToWomBundle: WomBundleMaker[FileElement] = new WomBundleMaker[FileElement] {

    override def toWomBundle(a: FileElement, importResolvers: List[String => Future[Checked[WomBundle]]]): Checked[WomBundle] = {

      val importsValidation: ErrorOr[Vector[ImportElement]] = if (a.imports.isEmpty) Vector.empty.valid else "FileElement to WOM conversion of imports not yet implemented.".invalidNel
      val tasksValidation: ErrorOr[Vector[TaskDefinitionElement]] = if (a.imports.isEmpty) Vector.empty.valid else "FileElement to WOM conversion of tasks not yet implemented.".invalidNel

      // TODO: Handle imports:
      val imports: Set[WomBundle] = Set.empty

      val structsValidation: ErrorOr[Map[String, WomType]] = StructEvaluation.convert(StructEvaluationInputs(a.structs, imports.flatMap(_.typeAliases).toMap))

      def toWorkflowInner(imports: Vector[ImportElement], tasks: Vector[TaskDefinitionElement], structs: Map[String, WomType]): ErrorOr[WomBundle] = {
        implicit val workflowConverter: CheckedAtoB[WorkflowDefinitionConvertInputs, WorkflowDefinition] = workflowDefinitionElementToWomWorkflowDefinition

        val tasks: Set[TaskDefinition] = Set.empty

        val workflowsValidation: ErrorOr[Vector[WorkflowDefinition]] = {
          a.workflows.toVector.traverse[ErrorOr, WorkflowDefinition] { workflowDefinition =>
            val convertInputs = WorkflowDefinitionConvertInputs(workflowDefinition, structs)
            workflowConverter.run(convertInputs).toValidated
          }
        }

        workflowsValidation map { workflows =>
          WomBundle(tasks ++ workflows, structs)
        }
      }

      ((importsValidation, tasksValidation, structsValidation) flatMapN { toWorkflowInner }).toEither
    }
  }

  def convert(a: FileElementAndImportResolvers): Checked[WomBundle] = a.fileElement.toWomBundle(a.importResolvers)
}

final case class FileElementAndImportResolvers(fileElement: FileElement, importResolvers: List[String => Future[Checked[WomBundle]]])