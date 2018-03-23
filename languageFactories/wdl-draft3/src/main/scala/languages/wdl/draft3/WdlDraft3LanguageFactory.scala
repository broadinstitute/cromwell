package languages.wdl.draft3

// TODO 2.11: Yup, cats.syntax.either._ again:
import java.nio.file.Paths

import better.files.File
import cats.syntax.either._
import cats.syntax.validated._
import cats.instances.either._
import cats.data.EitherT.fromEither
import cats.effect.IO
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import common.validation.Parse.Parse
import cromwell.core._
import cromwell.core.path.Path
import cromwell.languages.LanguageFactory.ImportResolver
import cromwell.languages.util.LanguageFactoryUtil
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import languages.wdl.draft3.WdlDraft3LanguageFactory._
import wdl.draft3.transforms.ast2wdlom._
import wdl.draft3.transforms.parsing._
import wdl.draft3.transforms.wdlom2wom._
import wdl.draft3.transforms.wdlom2wom.WomBundleToWomExecutable._
import wom.core.{WorkflowJson, WorkflowOptionsJson, WorkflowSource}
import wom.executable.WomBundle
import wom.expression.IoFunctionSet
import wom.transforms.WomExecutableMaker.ops._

import scala.util.Try

class WdlDraft3LanguageFactory() extends LanguageFactory {
  override def validateNamespace(source: WorkflowSourceFilesCollection,
                                    workflowOptions: WorkflowOptions,
                                    importLocalFilesystem: Boolean,
                                    workflowIdForLogging: WorkflowId,
                                    ioFunctions: IoFunctionSet): Parse[ValidatedWomNamespace] = {

    val factories: List[LanguageFactory] = List(this)
    val importResolvers: List[ImportResolver] = source.importsZipFileOption.map(zippedImportsResolver).toList

    val errorOr: Checked[ValidatedWomNamespace] = for {
      bundle <- getWomBundle(source.workflowSource, source.workflowOptionsJson, importResolvers, factories)
      executable <- createExecutable(bundle, source.inputsJson, ioFunctions)
    } yield executable

    fromEither[IO](errorOr)
  }

  override def getWomBundle(workflowSource: WorkflowSource, workflowOptionsJson: WorkflowOptionsJson, importResolvers: List[ImportResolver], languageFactories: List[LanguageFactory]): Checked[WomBundle] = {
    val converter: CheckedAtoB[FileStringParserInput, WomBundle] = stringToAst andThen astToFileElement.map(FileElementToWomBundleInputs(_, workflowOptionsJson, importResolvers, languageFactories)) andThen fileElementToWomBundle
    converter.run(FileStringParserInput(workflowSource, "input.wdl"))
  }

  override def createExecutable(womBundle: WomBundle, inputsJson: WorkflowJson, ioFunctions: IoFunctionSet): Checked[ValidatedWomNamespace] = {

    womBundle.toWomExecutable(Option(inputsJson), ioFunctions) flatMap { executable =>
      LanguageFactoryUtil.validateWomNamespace(executable)
    }
  }
}

object WdlDraft3LanguageFactory {
  def zippedImportsResolver(zippedImports: Array[Byte]): ImportResolver = {
    directoryResolver(LanguageFactoryUtil.validateImportsDirectory(zippedImports))
  }

  def directoryResolver(directoryValidation: ErrorOr[Path]): ImportResolver = CheckedAtoB.fromErrorOr { path =>
    directoryValidation flatMap { directory =>
      Try(Paths.get(directory.resolve(path).toFile.getCanonicalPath)).toErrorOr flatMap { absolutePathToFile =>
        val absolutePathToImports = Paths.get(directory.toJava.getCanonicalPath)
        if (absolutePathToFile.startsWith(absolutePathToImports)) {
          val file = File(absolutePathToFile)
          if (file.exists) {
            File(absolutePathToFile).contentAsString.validNel
          } else {
            s"Import file not found: $path".invalidNel
          }
        } else {
          s"$path is not a valid import".invalidNel
        }
      }
    }
  }
}
