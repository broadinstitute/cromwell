package wdl.draft3.transforms.wdlom2wom

import cats.instances.vector._
import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.traverse._
import cats.syntax.validated._
import cats.instances.either._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import common.validation.Checked._
import cromwell.languages.LanguageFactory
import cromwell.languages.LanguageFactory.ImportResolver
import wdl.draft3.transforms.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition.WorkflowDefinitionConvertInputs
import wdl.model.draft3.elements._
import wom.callable.{TaskDefinition, WorkflowDefinition}
import wom.executable.WomBundle
import wom.transforms.WomBundleMaker
import wom.transforms.WomBundleMaker.ops._
import wdl.draft3.transforms.wdlom2wom.StructEvaluation.StructEvaluationInputs
import wom.core.{WorkflowOptionsJson, WorkflowSource}
import wom.types.WomType

object FileElementToWomBundle {

  implicit val fileElementToWomBundle: WomBundleMaker[FileElementToWomBundleInputs] = new WomBundleMaker[FileElementToWomBundleInputs] {

    override def toWomBundle(a: FileElementToWomBundleInputs): Checked[WomBundle] = {

      val tasksValidation: ErrorOr[Vector[TaskDefinitionElement]] = a.fileElement.tasks.toVector.validNel
      val importsValidation: ErrorOr[Vector[WomBundle]] = a.fileElement.imports.toVector.traverse[ErrorOr, WomBundle] { importWomBundle(_, a.workflowOptionsJson, a.importResolvers, a.languageFactories) }

      def toWorkflowInner(imports: Vector[WomBundle], tasks: Vector[TaskDefinitionElement], structs: Map[String, WomType]): ErrorOr[WomBundle] = {
        val workflowConverter: CheckedAtoB[WorkflowDefinitionConvertInputs, WorkflowDefinition] = workflowDefinitionElementToWomWorkflowDefinition
        val taskConverter: CheckedAtoB[TaskDefinitionElement, TaskDefinition] = taskDefinitionElementToWomTaskDefinition

        val allStructs = structs ++ imports.flatMap(_.typeAliases)

        val taskDefs: ErrorOr[Set[TaskDefinition]] = {
          tasks.traverse[ErrorOr, TaskDefinition] { taskDefinition =>
            taskConverter.run(taskDefinition).toValidated
          }.map(_.toSet)
        }

        val workflowsValidation: ErrorOr[Vector[WorkflowDefinition]] = {
          a.fileElement.workflows.toVector.traverse[ErrorOr, WorkflowDefinition] { workflowDefinition =>
            val convertInputs = WorkflowDefinitionConvertInputs(workflowDefinition, allStructs)
            workflowConverter.run(convertInputs).toValidated
          }
        }

        (workflowsValidation, taskDefs) mapN { (workflows, tasks) =>
          WomBundle(tasks ++ workflows, allStructs)
        }
      }

      (importsValidation flatMap { imports =>
        val structsValidation: ErrorOr[Map[String, WomType]] = StructEvaluation.convert(StructEvaluationInputs(a.fileElement.structs, imports.flatMap(_.typeAliases).toMap))
        (tasksValidation, structsValidation) flatMapN { (tasks, structs) => toWorkflowInner(imports, tasks, structs) }
      }).toEither
    }
  }

  def convert(a: FileElementToWomBundleInputs): Checked[WomBundle] = a.toWomBundle

  private def importWomBundle(importElement: ImportElement,
                              optionsJson: WorkflowOptionsJson,
                              importResolvers: List[ImportResolver],
                              languageFactories: List[LanguageFactory]): ErrorOr[WomBundle] = {
    val compoundImportResolver: CheckedAtoB[String, WorkflowSource] = CheckedAtoB.firstSuccess(importResolvers, s"resolve import '${importElement.importUrl}'")

    val languageFactoryKleislis: List[CheckedAtoB[WorkflowSource, WomBundle]] = languageFactories map { factory =>
      CheckedAtoB.fromCheck { source: WorkflowSource =>
        factory.getWomBundle(source, optionsJson, importResolvers, languageFactories)
      }
    }
    val compoundLanguageFactory: CheckedAtoB[WorkflowSource, WomBundle] = CheckedAtoB.firstSuccess(languageFactoryKleislis, s"convert imported '${importElement.importUrl}' to WOM")

    val overallConversion = compoundImportResolver andThen compoundLanguageFactory

    (overallConversion.run(importElement.importUrl) flatMap { respectImportRenames(_, importElement.structRenames) }).toValidated
  }

  private def respectImportRenames(womBundle: WomBundle, importAliases: Map[String, String]): Checked[WomBundle] = {
    val importedStructs = womBundle.typeAliases
    val unexpectedAliases = importAliases.keySet.diff(womBundle.typeAliases.keySet)
    if (unexpectedAliases.isEmpty) {
      val newStructs = importedStructs map {
        case (key, value) if importAliases.contains(key) => importAliases(key) -> value
        case (otherKey, otherValue) => otherKey -> otherValue
      }
      womBundle.copy(typeAliases = newStructs).validNelCheck
    } else {
      s"Cannot import and rename: [${unexpectedAliases.mkString(", ")}] because the set of imported structs was: [${importedStructs.keySet.mkString(", ")}]".invalidNelCheck
    }
  }
}

final case class FileElementToWomBundleInputs(fileElement: FileElement,
                                              workflowOptionsJson: WorkflowOptionsJson,
                                              importResolvers: List[ImportResolver],
                                              languageFactories: List[LanguageFactory])
