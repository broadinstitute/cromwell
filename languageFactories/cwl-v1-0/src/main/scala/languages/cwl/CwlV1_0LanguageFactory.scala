package languages.cwl

import better.files.File
import cats.data.{EitherT, NonEmptyList}
import cats.data.EitherT.fromEither
import cats.effect.IO
import com.typesafe.config.Config
import common.Checked
import common.validation.Checked._
import common.validation.Parse.Parse
import cromwell.core.{WorkflowId, WorkflowOptions, WorkflowSourceFilesCollection}
import cromwell.languages.util.ImportResolver.ImportResolver
import cromwell.languages.util.LanguageFactoryUtil
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import cwl.{Cwl, CwlDecoder}
import wom.core.{WorkflowJson, WorkflowOptionsJson, WorkflowSource}
import wom.executable.{Executable, WomBundle}
import wom.expression.{IoFunctionSet, NoIoFunctionSet}

class CwlV1_0LanguageFactory(override val config: Config) extends LanguageFactory {

  override val languageName: String = "CWL"
  override val languageVersionName: String = "v1.0"

  override def validateNamespace(source: WorkflowSourceFilesCollection,
                                 workflowSource: WorkflowSource,
                                 workflowOptions: WorkflowOptions,
                                 importLocalFilesystem: Boolean,
                                 workflowIdForLogging: WorkflowId,
                                 ioFunctions: IoFunctionSet,
                                 importResolvers: List[ImportResolver]): Parse[ValidatedWomNamespace] = {
    import cwl.AcceptAllRequirements
    for {
      _ <- fromEither[IO](enabledCheck)
      cwl <- CwlDecoder.decodeCwlString(
        workflowSource,
        source.importsZipFileOption.map(File.newTemporaryFile().appendByteArray(_)),
        source.workflowRoot,
        "cwl_temp_file_" + workflowIdForLogging.toString
      )
      executable <- fromEither[IO](cwl.womExecutable(AcceptAllRequirements, Option(source.inputsJson), ioFunctions, strictValidation))
      validatedWomNamespace <- fromEither[IO](LanguageFactoryUtil.validateWomNamespace(executable, ioFunctions))
    } yield validatedWomNamespace
  }

  override def getWomBundle(workflowSource: WorkflowSource,
                            workflowOptionsJson: WorkflowOptionsJson,
                            importResolvers: List[ImportResolver],
                            languageFactories: List[LanguageFactory]): Checked[WomBundle] = {
    val cwlParse: Parse[Cwl] = CwlDecoder.decodeCwlString(workflowSource)

    import cwl.AcceptAllRequirements

    val executableIO: EitherT[IO, NonEmptyList[String], Executable] = for {
      cwl <- cwlParse
      executable <- fromEither[IO](cwl.womExecutable(AcceptAllRequirements, None, NoIoFunctionSet, strictValidation))
    } yield executable

    executableIO.value.unsafeRunSync() match {
      case Right(value) => WomBundle(Option(value.entryPoint), Map.empty, Map.empty).validNelCheck
      case Left(errors) =>
        val formattedErrors = errors.toList.mkString(System.lineSeparator(), System.lineSeparator(), System.lineSeparator())
        throw new Exception(formattedErrors)
    }

  }

  // wom.executable.WomBundle.toExecutableCallable
  // wom.executable.Executable
  override def createExecutable(womBundle: WomBundle, inputs: WorkflowJson, ioFunctions: IoFunctionSet): Checked[ValidatedWomNamespace] = {
    val x = womBundle.toExecutableCallable

    ValidatedWomNamespace(Executable(x.right.get, Map.empty), Map.empty, Map.empty).validNelCheck
  }

  override def looksParsable(content: String): Boolean = content.lines.exists { l =>
    val trimmed = l.trim.stripSuffix(",")
    trimmed == """"cwlVersion": "v1.0"""" || trimmed == "cwlVersion: v1.0"
  }
}
