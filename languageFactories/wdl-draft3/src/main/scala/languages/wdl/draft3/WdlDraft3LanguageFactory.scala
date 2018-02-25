package languages.wdl.draft3

import cats.data.EitherT.fromEither
import cats.effect.IO
import common.Checked
import common.validation.Parse.Parse
import cromwell.core._
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import wdl.draft3.transforms.ast2wdlom._
import wdl.draft3.transforms.parsing.{FileStringParserInput, StringParser}
import wdl.draft3.transforms.wdlom2wom._
import wom.transforms.WomExecutableMaker.ExecutableMakerInputs

class WdlDraft3LanguageFactory() extends LanguageFactory {
  override def validateNamespace(source: WorkflowSourceFilesCollection,
                                    workflowOptions: WorkflowOptions,
                                    importLocalFilesystem: Boolean,
                                    workflowIdForLogging: WorkflowId): Parse[ValidatedWomNamespace] = {


    val errorOr: Checked[ValidatedWomNamespace] = for {
      ast <- StringParser.convert(FileStringParserInput(source.workflowSource, "top-level.wdl"))
      fileElement <- astToFileElement(ast)
      executable <- fileElementToWomExecutable(ExecutableMakerInputs(fileElement, List.empty, Option(source.inputsJson)))

    } yield ValidatedWomNamespace(executable, Map.empty, Map.empty)

    fromEither[IO](errorOr)
  }

}
