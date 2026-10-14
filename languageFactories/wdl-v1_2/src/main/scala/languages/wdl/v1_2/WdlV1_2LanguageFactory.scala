package languages.wdl.v1_2

import java.util.concurrent.Callable

import cats.syntax.either._
import cats.instances.either._
import cats.data.EitherT.fromEither
import cats.effect.IO
import com.typesafe.config.Config
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import common.validation.IOChecked.IOChecked
import cromwell.core._
import cromwell.languages.util.ImportResolver._
import cromwell.languages.util.{LanguageFactoryUtil, ParserCache}
import cromwell.languages.util.ParserCache.ParserCacheInputs
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import wdl.transforms.base.wdlom2wom.WomBundleToWomExecutable._
import wdl.transforms.base.wdlom2wom._
import wdl.transforms.v1_2.ast2wdlom._
import wdl.transforms.v1_2.parsing._
import wdl.transforms.v1_2.wdlom2wom._
import wom.ResolvedImportRecord
import wom.core.{WorkflowJson, WorkflowOptionsJson, WorkflowSource}
import wom.executable.WomBundle
import wom.expression.IoFunctionSet
import wom.transforms.WomExecutableMaker.ops._

/**
 * Language factory for WDL 1.2 (https://github.com/openwdl/wdl/releases/tag/v1.2.1).
 *
 * Supported backends: GCP Batch, AWS Batch, TES and SFS (Local/HPC).
 * Note: WDL 1.2 is not supported yet on the PAPI v2 (Google Pipelines API / Life Sciences) backend. Use GCP Batch.
 *
 * WDL 1.2 documents can only import other WDL 1.2 documents, as the spec requires, so `factories` only contains this
 * factory.
 */
class WdlV1_2LanguageFactory(override val config: Config) extends LanguageFactory with ParserCache[WomBundle] {

  override val languageName: String = "WDL"
  override val languageVersionName: String = "1.2"

  override def validateNamespace(source: WorkflowSourceFilesCollection,
                                 workflowSource: WorkflowSource,
                                 workflowOptions: WorkflowOptions,
                                 importLocalFilesystem: Boolean,
                                 workflowIdForLogging: WorkflowId,
                                 ioFunctions: IoFunctionSet,
                                 importResolvers: List[ImportResolver]
  ): IOChecked[ValidatedWomNamespace] = {

    val factories: List[LanguageFactory] = List(this)

    val checked: Checked[ValidatedWomNamespace] = for {
      _ <- enabledCheck
      bundle <- getWomBundle(workflowSource,
                             workflowSourceOrigin = None,
                             source.workflowOptions.asPrettyJson,
                             importResolvers,
                             factories
      )
      executable <- createExecutable(bundle, source.inputsJson, ioFunctions)
    } yield executable

    fromEither[IO](checked)

  }

  override def getWomBundle(workflowSource: WorkflowSource,
                            workflowSourceOrigin: Option[ResolvedImportRecord],
                            workflowOptionsJson: WorkflowOptionsJson,
                            importResolvers: List[ImportResolver],
                            languageFactories: List[LanguageFactory],
                            convertNestedScatterToSubworkflow: Boolean = true
  ): Checked[WomBundle] = {

    val converter: CheckedAtoB[FileStringParserInput, WomBundle] =
      stringToAst andThen wrapAst andThen astToFileElement.map(
        FileElementToWomBundleInputs(
          _,
          workflowOptionsJson,
          convertNestedScatterToSubworkflow,
          importResolvers,
          languageFactories,
          workflowDefinitionElementToWomWorkflowDefinition,
          taskDefinitionElementToWomTaskDefinition
        )
      ) andThen fileElementToWomBundle

    lazy val validationCallable = new Callable[ErrorOr[WomBundle]] {
      def call: ErrorOr[WomBundle] = converter
        .run(FileStringParserInput(workflowSource, workflowSourceOrigin.map(_.importPath).getOrElse("input.wdl")))
        .map(b => b.copyResolvedImportRecord(b, workflowSourceOrigin))
        .toValidated
    }

    lazy val parserCacheInputs =
      ParserCacheInputs(Option(workflowSource), workflowSourceOrigin.map(_.importPath), None, importResolvers)

    for {
      _ <- enabledCheck
      womBundle <- retrieveOrCalculate(parserCacheInputs, validationCallable).toEither
    } yield womBundle
  }

  // TODO: WDL 1.2 recommends an "extended" JSON format for File and Directory inputs and outputs (objects with
  //  `location`, `basename` and a `listing` of a directory's contents), which helps reproducibility and call caching.
  //  Only the plain string form is accepted today. Suggestion: support the extended input format when coercing input
  //  JSON to File/Directory values, and decide whether the outputs/metadata API should return the extended format for
  //  WDL 1.2 workflows (it changes the shape of the outputs JSON).
  // TODO: WDL 1.2.1 requires File and Directory values to be canonicalized and to exist when they are created, and a
  //  File must not refer to a directory (or vice versa). Cromwell checks existence when files are localized. Checking
  //  cloud paths when every declaration is evaluated costs extra API calls. Suggestion: validate workflow inputs here,
  //  and leave other declarations to localization.
  override def createExecutable(womBundle: WomBundle,
                                inputsJson: WorkflowJson,
                                ioFunctions: IoFunctionSet
  ): Checked[ValidatedWomNamespace] =
    for {
      _ <- enabledCheck
      executable <- womBundle.toWomExecutable(Option(inputsJson), ioFunctions, strictValidation)
      validated <- LanguageFactoryUtil.validateWomNamespace(executable, ioFunctions)
    } yield validated

  override def looksParsable(content: String): Boolean =
    LanguageFactoryUtil.simpleLooksParseable(List("version 1.2"), List("#"))(content)
}
