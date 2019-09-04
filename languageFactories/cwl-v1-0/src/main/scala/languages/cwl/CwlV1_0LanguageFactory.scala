package languages.cwl

import better.files.File
import cats.data.EitherT.fromEither
import cats.data.NonEmptyList
import cats.effect.IO
import com.typesafe.config.Config
import common.Checked
import common.validation.Checked._
import common.validation.IOChecked.IOChecked
import cromwell.core.{WorkflowId, WorkflowOptions, WorkflowSourceFilesCollection}
import cromwell.languages.util.ImportResolver.ImportResolver
import cromwell.languages.util.LanguageFactoryUtil
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import cwl.preprocessor.CwlReference
import cwl.{Cwl, CwlDecoder}
import wom.ResolvedImportRecord
import wom.core.{WorkflowJson, WorkflowOptionsJson, WorkflowSource}
import wom.executable.WomBundle
import wom.expression.IoFunctionSet

class CwlV1_0LanguageFactory(override val config: Config) extends LanguageFactory {

  override val languageName: String = "CWL"
  override val languageVersionName: String = "v1.0"

  override def validateNamespace(source: WorkflowSourceFilesCollection,
                                 workflowSource: WorkflowSource,
                                 workflowOptions: WorkflowOptions,
                                 importLocalFilesystem: Boolean,
                                 workflowIdForLogging: WorkflowId,
                                 ioFunctions: IoFunctionSet,
                                 importResolvers: List[ImportResolver]): IOChecked[ValidatedWomNamespace] = {

    def parse(): IOChecked[Cwl] = source.workflowUrl match {
      case Some(url) => for {
        reference <- fromEither[IO](CwlReference.fromString(url).map(Right(_)).getOrElse(Left(NonEmptyList.one(s"Invalid workflow reference: $url")))): IOChecked[CwlReference]
        parsed <- CwlDecoder.decodeCwlReference(reference.changePointer(source.workflowRoot))
      } yield parsed
      case None =>
        CwlDecoder.decodeCwlString(
          workflowSource,
          source.importsZipFileOption.map(File.newTemporaryFile().appendByteArray(_)),
          source.workflowRoot,
          "cwl_temp_file_" + workflowIdForLogging.toString
        )
    }

    import cwl.AcceptAllRequirements
    for {
      _ <- fromEither[IO](enabledCheck)
      cwl <- parse()
      executable <- fromEither[IO](cwl.womExecutable(AcceptAllRequirements, Option(source.inputsJson), ioFunctions, strictValidation))
      validatedWomNamespace <- fromEither[IO](LanguageFactoryUtil.validateWomNamespace(executable, ioFunctions))
    } yield validatedWomNamespace
  }

  override def getWomBundle(workflowSource: WorkflowSource,
                            workflowSourceOrigin: Option[ResolvedImportRecord],
                            workflowOptionsJson: WorkflowOptionsJson,
                            importResolvers: List[ImportResolver],
                            languageFactories: List[LanguageFactory],
                            convertNestedScatterToSubworkflow : Boolean = true): Checked[WomBundle] =
    enabledCheck flatMap { _ => "No getWomBundle method implemented in CWL v1".invalidNelCheck }

  override def createExecutable(womBundle: WomBundle, inputs: WorkflowJson, ioFunctions: IoFunctionSet): Checked[ValidatedWomNamespace] =
    enabledCheck flatMap { _ => "No createExecutable method implemented in CWL v1".invalidNelCheck }

  override def looksParsable(content: String): Boolean = content.linesIterator.exists { l =>
    val trimmed = l.trim.stripSuffix(",")
    trimmed == """"cwlVersion": "v1.0"""" || trimmed == "cwlVersion: v1.0"
  }
}
