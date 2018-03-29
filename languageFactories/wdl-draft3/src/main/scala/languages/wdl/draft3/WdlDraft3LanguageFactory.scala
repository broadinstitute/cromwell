package languages.wdl.draft3

// TODO 2.11: Yup, cats.syntax.either._ again:
import cats.syntax.either._
import cats.instances.either._
import cats.data.EitherT.fromEither
import cats.effect.IO
import common.Checked
import common.validation.Checked._
import common.transforms.CheckedAtoB
import common.validation.Parse.Parse
import cromwell.core._
import cromwell.languages.util.ImportResolver._
import cromwell.languages.util.LanguageFactoryUtil
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import wdl.draft3.transforms.ast2wdlom._
import wdl.draft3.transforms.parsing._
import wdl.draft3.transforms.wdlom2wom._
import wdl.draft3.transforms.wdlom2wom.WomBundleToWomExecutable._
import wom.core.{WorkflowJson, WorkflowOptionsJson, WorkflowSource}
import wom.executable.WomBundle
import wom.expression.IoFunctionSet
import wom.transforms.WomExecutableMaker.ops._

class WdlDraft3LanguageFactory(override val config: Map[String, Any]) extends LanguageFactory {

  override def validateNamespace(source: WorkflowSourceFilesCollection,
                                    workflowOptions: WorkflowOptions,
                                    importLocalFilesystem: Boolean,
                                    workflowIdForLogging: WorkflowId,
                                    ioFunctions: IoFunctionSet): Parse[ValidatedWomNamespace] = {

    val factories: List[LanguageFactory] = List(this)
    val localFilesystemResolvers = if (importLocalFilesystem) List(localFileResolver) else List.empty
    val importResolvers: List[ImportResolver] = source.importsZipFileOption.map(zippedImportsResolver).toList ++ localFilesystemResolvers

    val errorOr: Checked[ValidatedWomNamespace] = for {
      _ <- standardConfig.enabledCheck
      bundle <- getWomBundle(source.workflowSource, source.workflowOptionsJson, importResolvers, factories)
      executable <- createExecutable(bundle, source.inputsJson, ioFunctions)
    } yield executable

    fromEither[IO](errorOr)

  }

  override def getWomBundle(workflowSource: WorkflowSource, workflowOptionsJson: WorkflowOptionsJson, importResolvers: List[ImportResolver], languageFactories: List[LanguageFactory]): Checked[WomBundle] = {
    val checkEnabled: CheckedAtoB[FileStringParserInput, FileStringParserInput] = CheckedAtoB.fromCheck(x => standardConfig.enabledCheck map(_ => x))
    val converter: CheckedAtoB[FileStringParserInput, WomBundle] = checkEnabled andThen stringToAst andThen astToFileElement.map(FileElementToWomBundleInputs(_, workflowOptionsJson, importResolvers, languageFactories)) andThen fileElementToWomBundle
    converter.run(FileStringParserInput(workflowSource, "input.wdl"))
  }

  override def createExecutable(womBundle: WomBundle, inputsJson: WorkflowJson, ioFunctions: IoFunctionSet): Checked[ValidatedWomNamespace] = {
    for {
      _ <- standardConfig.enabledCheck
      executable <- womBundle.toWomExecutable(Option(inputsJson), ioFunctions, standardConfig.strictValidation)
      validated <- LanguageFactoryUtil.validateWomNamespace(executable)
    } yield validated
  }
}
