package languages.wdl.biscayne

import java.util.concurrent.Callable

import cats.syntax.either._
import cats.data.EitherT.fromEither
import cats.effect.IO
import cats.instances.either._
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
import wdl.transforms.biscayne.ast2wdlom._
import wdl.transforms.biscayne.parsing._
import wdl.transforms.biscayne.wdlom2wom._
import wom.ResolvedImportRecord
import wom.core.{WorkflowJson, WorkflowOptionsJson, WorkflowSource}
import wom.executable.WomBundle
import wom.expression.IoFunctionSet
import wom.transforms.WomExecutableMaker.ops._

class WdlBiscayneLanguageFactory(override val config: Config) extends LanguageFactory with ParserCache[WomBundle] {

  override val languageName: String = "WDL"
  override val languageVersionName: String = "Biscayne"

  override def validateNamespace(source: WorkflowSourceFilesCollection,
                                 workflowSource: WorkflowSource,
                                 workflowOptions: WorkflowOptions,
                                 importLocalFilesystem: Boolean,
                                 workflowIdForLogging: WorkflowId,
                                 ioFunctions: IoFunctionSet,
                                 importResolvers: List[ImportResolver]): IOChecked[ValidatedWomNamespace] = {

    val factories: List[LanguageFactory] = List(this)

    val checked: Checked[ValidatedWomNamespace] = for {
      _ <- enabledCheck
      bundle <- getWomBundle(workflowSource, workflowSourceOrigin = None, source.workflowOptions.asPrettyJson, importResolvers, factories)
      executable <- createExecutable(bundle, source.inputsJson, ioFunctions)
    } yield executable

    fromEither[IO](checked)

  }

  override def getWomBundle(workflowSource: WorkflowSource,
                            workflowSourceOrigin: Option[ResolvedImportRecord],
                            workflowOptionsJson: WorkflowOptionsJson,
                            importResolvers: List[ImportResolver],
                            languageFactories: List[LanguageFactory],
                            convertNestedScatterToSubworkflow : Boolean = true): Checked[WomBundle] = {

    val converter: CheckedAtoB[FileStringParserInput, WomBundle] = stringToAst andThen wrapAst andThen astToFileElement.map(FileElementToWomBundleInputs(_, workflowOptionsJson, convertNestedScatterToSubworkflow, importResolvers, languageFactories, workflowDefinitionElementToWomWorkflowDefinition, taskDefinitionElementToWomTaskDefinition)) andThen fileElementToWomBundle

    lazy val validationCallable = new Callable[ErrorOr[WomBundle]] {
      def call: ErrorOr[WomBundle] = converter
        .run(FileStringParserInput(workflowSource, workflowSourceOrigin.map(_.importPath).getOrElse("input.wdl")))
        .map(b => b.copyResolvedImportRecord(b, workflowSourceOrigin)).toValidated
    }

    lazy val parserCacheInputs = ParserCacheInputs(Option(workflowSource), workflowSourceOrigin.map(_.importPath), None, importResolvers)

    for {
      _ <- enabledCheck
      womBundle <- retrieveOrCalculate(parserCacheInputs, validationCallable).toEither
    } yield womBundle
  }

  override def createExecutable(womBundle: WomBundle, inputsJson: WorkflowJson, ioFunctions: IoFunctionSet): Checked[ValidatedWomNamespace] = {
    for {
      _ <- enabledCheck
      executable <- womBundle.toWomExecutable(Option(inputsJson), ioFunctions, strictValidation)
      validated <- LanguageFactoryUtil.validateWomNamespace(executable, ioFunctions)
    } yield validated
  }

  override def looksParsable(content: String): Boolean = LanguageFactoryUtil.simpleLooksParseable(List("version development"), List("#"))(content)
}
