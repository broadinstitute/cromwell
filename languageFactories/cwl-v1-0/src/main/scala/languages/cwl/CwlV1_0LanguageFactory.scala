package languages.cwl

import better.files.File
import cats.Monad
import cats.data.EitherT.fromEither
import cats.effect.IO
import common.Checked
import common.validation.Checked._
import common.validation.Parse.{Parse, errorOrParse, goParse, tryParse}
import cromwell.core.path.DefaultPathBuilder
import cromwell.core.{WorkflowId, WorkflowOptions, WorkflowSourceFilesCollection, WorkflowSourceFilesWithDependenciesZip}
import cromwell.languages.util.ImportResolver.ImportResolver
import cromwell.languages.util.LanguageFactoryUtil
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import cwl.CwlDecoder
import wom.core.{WorkflowJson, WorkflowOptionsJson, WorkflowSource}
import wom.executable.WomBundle
import wom.expression.IoFunctionSet

class CwlV1_0LanguageFactory(override val config: Map[String, Any]) extends LanguageFactory {

  override val languageName: String = "CWL"
  override val languageVersionName: String = "v1.0"

  override def validateNamespace(source: WorkflowSourceFilesCollection,
                                 workflowOptions: WorkflowOptions,
                                 importLocalFilesystem: Boolean,
                                 workflowIdForLogging: WorkflowId,
                                 ioFunctions: IoFunctionSet): Parse[ValidatedWomNamespace] = {
    // TODO WOM: CwlDecoder takes a file so write it to disk for now

    def writeCwlFileToNewTempDir(): Parse[File] = {
      goParse {
        val tempDir = File.newTemporaryDirectory(prefix = s"$workflowIdForLogging.temp.")
        val cwlFile: File = tempDir./(s"$workflowIdForLogging.cwl").write(source.workflowSource)
        cwlFile
      }
    }

    def unzipDependencies(cwlFile: File): Parse[Unit] = {
      source match {
        case wsfwdz: WorkflowSourceFilesWithDependenciesZip =>
          for {
            parent <- tryParse(DefaultPathBuilder.build(cwlFile.parent.pathAsString))
            _ <- errorOrParse(LanguageFactoryUtil.validateImportsDirectory(wsfwdz.importsZip, Option(parent)))
          } yield ()
        case _ => Monad[Parse].unit
      }
    }

    import cwl.AcceptAllRequirements
    for {
      _ <- fromEither[IO](standardConfig.enabledCheck)
      cwlFile <- writeCwlFileToNewTempDir()
      _ <- unzipDependencies(cwlFile)
      cwl <- CwlDecoder.decodeCwlFile(cwlFile, source.workflowRoot)
      executable <- fromEither[IO](cwl.womExecutable(AcceptAllRequirements, Option(source.inputsJson), ioFunctions, standardConfig.strictValidation))
      validatedWomNamespace <- fromEither[IO](LanguageFactoryUtil.validateWomNamespace(executable, ioFunctions))
      _ <- CwlDecoder.todoDeleteCwlFileParentDirectory(cwlFile.parent)
    } yield validatedWomNamespace
  }

  override def getWomBundle(workflowSource: WorkflowSource, workflowOptionsJson: WorkflowOptionsJson, importResolvers: List[ImportResolver], languageFactories: List[LanguageFactory]): Checked[WomBundle] =
    standardConfig.enabledCheck flatMap { _ => "No getWomBundle method implemented in CWL v1".invalidNelCheck }

  override def createExecutable(womBundle: WomBundle, inputs: WorkflowJson, ioFunctions: IoFunctionSet): Checked[ValidatedWomNamespace] =
    standardConfig.enabledCheck flatMap { _ => "No createExecutable method implemented in CWL v1".invalidNelCheck }

  override def looksParsable(content: String): Boolean = content.lines.exists { l =>
    val trimmed = l.trim.stripSuffix(",")
    trimmed == """"cwlVersion": "v1.0"""" || trimmed == "cwlVersion: v1.0"
  }
}
