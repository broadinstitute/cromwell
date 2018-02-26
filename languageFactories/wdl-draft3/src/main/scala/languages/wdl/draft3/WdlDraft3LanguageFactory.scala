package languages.wdl.draft3

// TODO 2.11: Yup, cats.syntax.either._ again:
import cats.syntax.either._
import cats.instances.either._
import cats.data.EitherT.fromEither
import cats.effect.IO
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.Parse.Parse
import cromwell.core._
import cromwell.languages.util.LanguageFactoryUtil
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import wdl.draft3.transforms.ast2wdlom._
import wdl.draft3.transforms.parsing._
import wdl.draft3.transforms.wdlom2wom._
import wdl.draft3.transforms.wdlom2wom.WomBundleToWomExecutable._
import wom.core.{WorkflowJson, WorkflowOptionsJson, WorkflowSource}
import wom.executable.WomBundle
import wom.transforms.WomExecutableMaker.ops._

import scala.concurrent.Future

class WdlDraft3LanguageFactory() extends LanguageFactory {
  override def validateNamespace(source: WorkflowSourceFilesCollection,
                                    workflowOptions: WorkflowOptions,
                                    importLocalFilesystem: Boolean,
                                    workflowIdForLogging: WorkflowId): Parse[ValidatedWomNamespace] = {


    val errorOr: Checked[ValidatedWomNamespace] = for {
      // TODO: Make a real list of import resolvers for draft 3:
      bundle <- getWomBundle(source.workflowSource, source.workflowOptionsJson, List.empty)
      executable <- createExecutable(bundle, source.inputsJson)
    } yield executable

    fromEither[IO](errorOr)
  }

  override def getWomBundle(workflowSource: WorkflowSource, workflowOptionsJson: WorkflowOptionsJson, importResolvers: List[String => Future[Checked[WomBundle]]]): Checked[WomBundle] = {
    val converter: CheckedAtoB[FileStringParserInput, WomBundle] = stringToAst andThen astToFileElement.map(FileElementAndImportResolvers(_, importResolvers)) andThen fileElementToWomBundle
    converter.run(FileStringParserInput(workflowSource, "input.wdl"))
  }

  override def createExecutable(womBundle: WomBundle, inputsJson: WorkflowJson): Checked[ValidatedWomNamespace] = {

    womBundle.toWomExecutable(Option(inputsJson)) flatMap { executable =>
      LanguageFactoryUtil.validateWomNamespace(executable)
    }
  }
}
