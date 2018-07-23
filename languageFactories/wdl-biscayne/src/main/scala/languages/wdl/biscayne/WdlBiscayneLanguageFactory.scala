package languages.wdl.biscayne

import java.nio.file.Paths

import cats.instances.either._
import cats.data.EitherT.fromEither
import cats.effect.IO
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.Parse.Parse
import common.validation.Checked._
import cromwell.core._
import cromwell.core.path.DefaultPath
import cromwell.languages.util.ImportResolver._
import cromwell.languages.util.LanguageFactoryUtil
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import wdl.transforms.biscayne.ast2wdlom._
import wdl.transforms.biscayne.parsing._
import wdl.transforms.base.wdlom2wom._
import wdl.transforms.base.wdlom2wom.WomBundleToWomExecutable._
import wom.core.{WorkflowJson, WorkflowOptionsJson, WorkflowSource}
import wom.executable.WomBundle
import wom.expression.IoFunctionSet
import wom.transforms.WomExecutableMaker.ops._
import wdl.transforms.biscayne.wdlom2wom._

class WdlBiscayneLanguageFactory(override val config: Map[String, Any]) extends LanguageFactory {

  override val languageName: String = "WDL"
  override val languageVersionName: String = "Biscayne"

  override def validateNamespace(source: WorkflowSourceFilesCollection,
                                    workflowOptions: WorkflowOptions,
                                    importLocalFilesystem: Boolean,
                                    workflowIdForLogging: WorkflowId,
                                    ioFunctions: IoFunctionSet): Parse[ValidatedWomNamespace] = {

    val factories: List[LanguageFactory] = List(this)
    val localFilesystemResolvers = if (importLocalFilesystem) List(DirectoryResolver(DefaultPath(Paths.get("/")))) else List.empty

    val zippedResolverCheck: Checked[Option[ImportResolver]] = source.importsZipFileOption match {
      case None => None.validNelCheck
      case Some(zipContent) => zippedImportResolver(zipContent).toEither.map(Option.apply)
    }

    val checked: Checked[ValidatedWomNamespace] = for {
      _ <- standardConfig.enabledCheck
      zippedImportResolver <- zippedResolverCheck
      importResolvers = zippedImportResolver.toList ++ localFilesystemResolvers :+ HttpResolver(None, Map.empty)
      bundle <- getWomBundle(source.workflowSource, source.workflowOptionsJson, importResolvers, factories)
      executable <- createExecutable(bundle, source.inputsJson, ioFunctions)
    } yield executable

    fromEither[IO](checked)

  }

  override def getWomBundle(workflowSource: WorkflowSource, workflowOptionsJson: WorkflowOptionsJson, importResolvers: List[ImportResolver], languageFactories: List[LanguageFactory]): Checked[WomBundle] = {
    val checkEnabled: CheckedAtoB[FileStringParserInput, FileStringParserInput] = CheckedAtoB.fromCheck(x => enabledCheck map(_ => x))
    val converter: CheckedAtoB[FileStringParserInput, WomBundle] = checkEnabled andThen stringToAst andThen wrapAst andThen astToFileElement.map(FileElementToWomBundleInputs(_, workflowOptionsJson, importResolvers, languageFactories, workflowDefinitionElementToWomWorkflowDefinition, taskDefinitionElementToWomTaskDefinition)) andThen fileElementToWomBundle
    converter.run(FileStringParserInput(workflowSource, "input.wdl"))
  }

  override def createExecutable(womBundle: WomBundle, inputsJson: WorkflowJson, ioFunctions: IoFunctionSet): Checked[ValidatedWomNamespace] = {
    for {
      _ <- enabledCheck
      executable <- womBundle.toWomExecutable(Option(inputsJson), ioFunctions, standardConfig.strictValidation)
      validated <- LanguageFactoryUtil.validateWomNamespace(executable, ioFunctions)
    } yield validated
  }

  override def looksParsable(content: String): Boolean = LanguageFactoryUtil.simpleLooksParseable(List("version 1.1", "version biscayne"), List("#"))(content)
}
