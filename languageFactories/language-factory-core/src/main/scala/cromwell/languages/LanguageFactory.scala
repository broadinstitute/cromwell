package cromwell.languages

import cats.data.NonEmptyList
import cats.data.Validated.Invalid
import cats.syntax.validated._
import cats.syntax.traverse._
import cats.syntax.functor._
import cats.instances.list._
import com.typesafe.config.Config
import common.Checked
import common.validation.Checked._
import common.validation.IOChecked.IOChecked
import cromwell.core.{WorkflowId, WorkflowOptions, WorkflowSourceFilesCollection}
import cromwell.languages.util.ImportResolver.ImportResolver
import wom.ResolvedImportRecord
import wom.callable.TaskDefinition
import wom.core._
import wom.executable.WomBundle
import wom.expression.IoFunctionSet
import wom.runtime.WomOutputRuntimeExtractor
import wom.types.WomSingleFileType

trait LanguageFactory {

  def languageName: String
  def languageVersionName: String

  // Passed in by the constructor:
  def config: Config

  import net.ceedubs.ficus.Ficus._

  lazy val enabled = !config.as[Option[Boolean]]("enabled").contains(false)
  lazy val enabledCheck: Checked[Unit] = if (enabled) ().validNelCheck else
    s"The language factory for $languageName ($languageVersionName) is not currently enabled in this Cromwell".invalidNelCheck


  lazy val strictValidation: Boolean = !config.as[Option[Boolean]]("strict-validation").contains(false)

  lazy val womOutputRuntimeExtractor: Checked[Option[WomOutputRuntimeExtractor]] = config.getAs[Config]("output-runtime-extractor") match {
    case Some(c) => WomOutputRuntimeExtractor.fromConfig(c).map(Option.apply).toEither
    case _ => None.validNelCheck
  }

  val allowOutputsAsFunctionsOfFileInputs = config.as[Option[Boolean]]("allow-outputs-as-functions-of-file-inputs").getOrElse(true)

  final def getWomBundle(workflowSource: WorkflowSource,
                         workflowSourceOrigin: Option[ResolvedImportRecord],
                         workflowOptionsJson: WorkflowOptionsJson,
                         importResolvers: List[ImportResolver],
                         languageFactories: List[LanguageFactory],
                         convertNestedScatterToSubworkflow : Boolean = true): Checked[WomBundle] = {
    val bundle = getWomBundleInner(workflowSource, workflowSourceOrigin, workflowOptionsJson, importResolvers, languageFactories, convertNestedScatterToSubworkflow)

    if (allowOutputsAsFunctionsOfFileInputs) bundle else {
      def validateOutputsAreNotFunctionsOfFiles(womBundle: WomBundle): Checked[Unit] = {
        val tasksToCheck = womBundle.allCallables.collect {
          case (name, t: TaskDefinition) => name -> t
        }.toList

        val taskValidation = tasksToCheck.traverse { case (name, task) =>
          def isAFileInput(inputName: String) = task.inputs.exists(i => i.localName.value == inputName && i.womType == WomSingleFileType)

          val invalidTasksOutputs = task.outputs.collect {
            case o if o.expression.inputs.exists(isAFileInput) => s"Cannot evaluate task '$name''s output '${o.name} = ${o.expression.sourceString}': it relies on File inputs."
          }

          NonEmptyList.fromList(invalidTasksOutputs) match {
            case None => ().validNel
            case Some(errors) => Invalid(errors)
          }
        }

        taskValidation.void.toEither
      }

      for {
        b <- bundle
        _ <- validateOutputsAreNotFunctionsOfFiles(b)
      } yield b
    }


  }

  protected def getWomBundleInner(workflowSource: WorkflowSource,
                                  workflowSourceOrigin: Option[ResolvedImportRecord],
                                  workflowOptionsJson: WorkflowOptionsJson,
                                  importResolvers: List[ImportResolver],
                                  languageFactories: List[LanguageFactory],
                                  convertNestedScatterToSubworkflow : Boolean = true): Checked[WomBundle]

  def createExecutable(womBundle: WomBundle,
                       inputs: WorkflowJson,
                       ioFunctions: IoFunctionSet): Checked[ValidatedWomNamespace]

  def validateNamespace(source: WorkflowSourceFilesCollection,
                        workflowSource: WorkflowSource,
                        workflowOptions: WorkflowOptions,
                        importLocalFilesystem: Boolean,
                        workflowIdForLogging: WorkflowId,
                        ioFunctions: IoFunctionSet,
                        importResolvers: List[ImportResolver]): IOChecked[ValidatedWomNamespace]

  /**
    * In case no version is specified: does this language factory feel like it might be suitable for this file?
    * @param content The workflow description
    */
  def looksParsable(content: String): Boolean
}
