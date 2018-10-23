package languages.cwl

import better.files.File
import cats.data.{EitherT, NonEmptyList}
import wom.callable.Callable
import cats.data.EitherT.fromEither
import cats.effect.IO
import com.typesafe.config.Config
import common.Checked
import common.validation.Checked._
import common.validation.Parse.Parse
import cromwell.core.path.Path
import cromwell.core.{WorkflowId, WorkflowOptions, WorkflowSourceFilesCollection}
import cromwell.languages.util.ImportResolver.ImportResolver
import cromwell.languages.util.{ImportResolver, LanguageFactoryUtil}
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import cwl.{Cwl, CwlDecoder}
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
    // Use the resolvers to snag the contents of specific files relative to the root workflow
    // Problem: file name to snag is not available until after CWL decoding (catch-22)
    //    val x: CheckedAtoB[ImportResolutionRequest, ImportResolver.ResolvedImportBundle] = CheckedAtoB.firstSuccess(importResolvers.map(_.resolver), "a")
    //    val y: Checked[ImportResolver.ResolvedImportBundle] = x.run(ImportResolutionRequest("bogus.cwl", importResolvers))
    //    val z = y.getOrElse(???).source
    //    println(z)

    val rootDir: Path = importResolvers(2).asInstanceOf[ImportResolver.DirectoryResolver].directory

//    import better.files.File.apply
    val rootDirZipped = File(rootDir.zip().toFile.toPath)

    val cwlParseAttempt: Parse[Cwl] = CwlDecoder.decodeCwlString(workflowSource, Option(rootDirZipped))

    import cwl.AcceptAllRequirements

    val callableIO: EitherT[IO, NonEmptyList[String], Either[NonEmptyList[String], Callable]] = for {
      cwl <- cwlParseAttempt
    } yield {
      cwl match {
        case Cwl.Workflow(workflow) =>
          workflow.womDefinition(AcceptAllRequirements, Vector.empty)
        case Cwl.CommandLineTool(cliTool) =>
          cliTool.buildTaskDefinition(AcceptAllRequirements, Vector.empty)
        case Cwl.ExpressionTool(expTool) =>
          expTool.buildTaskDefinition(AcceptAllRequirements, Vector.empty)
      }
    }

    callableIO.value.unsafeRunSync().joinRight match {
      case Right(value: Callable) => WomBundle(Option(value), Map.empty, Map.empty).validNelCheck
      case Left(errors) => errors.toList.mkString(", ").invalidNelCheck
    }
  }

  override def createExecutable(womBundle: WomBundle, inputs: WorkflowJson, ioFunctions: IoFunctionSet): Checked[ValidatedWomNamespace] =
    enabledCheck flatMap { _ => "No createExecutable method implemented in CWL v1".invalidNelCheck }

  override def looksParsable(content: String): Boolean = content.lines.exists { l =>
    val trimmed = l.trim.stripSuffix(",")
    trimmed == """"cwlVersion": "v1.0"""" || trimmed == "cwlVersion: v1.0"
  }
}
