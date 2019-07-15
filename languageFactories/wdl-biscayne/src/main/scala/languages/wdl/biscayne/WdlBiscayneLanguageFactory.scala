package languages.wdl.biscayne

import cats.data.EitherT.fromEither
import cats.effect.IO
import cats.instances.either._
import com.typesafe.config.Config
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.IOChecked.IOChecked
import cromwell.core._
import cromwell.languages.util.ImportResolver._
import cromwell.languages.util.LanguageFactoryUtil
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import wdl.transforms.base.wdlom2wom.WomBundleToWomExecutable._
import wdl.transforms.base.wdlom2wom._
import wdl.transforms.biscayne.ast2wdlom._
import wdl.transforms.biscayne.parsing._
import wdl.transforms.biscayne.wdlom2wom._
import wom.core.{WorkflowJson, WorkflowOptionsJson, WorkflowSource}
import wom.executable.WomBundle
import wom.expression.IoFunctionSet
import wom.transforms.WomExecutableMaker.ops._

class WdlBiscayneLanguageFactory(override val config: Config) extends LanguageFactory {

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
      bundle <- getWomBundle(workflowSource, source.workflowOptions.asPrettyJson, importResolvers, factories)
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
      executable <- womBundle.toWomExecutable(Option(inputsJson), ioFunctions, strictValidation)
      validated <- LanguageFactoryUtil.validateWomNamespace(executable, ioFunctions)
    } yield validated
  }

  override def looksParsable(content: String): Boolean = LanguageFactoryUtil.simpleLooksParseable(List("version development"), List("#"))(content)
}
