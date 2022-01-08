package cromwell.languages.util

import cats.data.NonEmptyList
import cats.syntax.validated._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import cromwell.core.CromwellGraphNode.CromwellEnhancedOutputPort
import cromwell.core.{WorkflowId, WorkflowSourceFilesCollection}
import cromwell.core.path.BetterFileMethods.OpenOptions
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.languages.config.CromwellLanguages
import cromwell.languages.util.ImportResolver.{ImportResolutionRequest, ImportResolver, ResolvedImportBundle}
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import wom.core.{WorkflowSource, WorkflowUrl}
import wom.executable.Executable
import wom.executable.Executable.ResolvedExecutableInputs
import wom.expression.IoFunctionSet
import wom.graph.GraphNodePort.OutputPort
import wom.values.{WomSingleFile, WomValue}

import scala.util.{Failure, Success, Try}

object LanguageFactoryUtil {

  /**
    * Unzip the imports.zip and validate that it was good.
    * @param zipContents the zip contents
    * @return where the imports were unzipped to
    */
  def createImportsDirectory(zipContents: Array[Byte], workflowId: WorkflowId): ErrorOr[Path] = {

    def makeZipFile: Try[Path] = Try {
      DefaultPathBuilder.createTempFile(s"imports_workflow_${workflowId}_", ".zip").writeByteArray(zipContents)(OpenOptions.default)
    }

    def unZipFile(f: Path) = Try(f.unzip)

    val importsFile = for {
      zipFile <- makeZipFile
      unzipped <- unZipFile(zipFile)
      _ <- Try(zipFile.delete(swallowIOExceptions = true))
    } yield unzipped

    importsFile match {
      case Success(unzippedDirectory: Path) => unzippedDirectory.validNel
      case Failure(t) => t.getMessage.invalidNel
    }
  }

  def validateWomNamespace(womExecutable: Executable, ioFunctions: IoFunctionSet): Checked[ValidatedWomNamespace] = for {
    evaluatedInputs <- validateExecutableInputs(womExecutable.resolvedExecutableInputs, ioFunctions).toEither
    validatedWomNamespace = ValidatedWomNamespace(womExecutable, evaluatedInputs, Map.empty)
    _ <- validateWdlFiles(validatedWomNamespace.womValueInputs)
  } yield validatedWomNamespace

  /*
    * At this point input values are either a WomValue (if it was provided through the input file)
    * or a WomExpression (if we fell back to the default).
    * We assume that default expressions do NOT reference any "variables" (other inputs, call outputs ...)
    * Should this assumption prove not sufficient InstantiatedExpressions or ExpressionNodes would need to be provided
    * instead so that they can be evaluated JIT.
    * Note that the ioFunctions use engine level pathBuilders. This means that their credentials come from the engine section
    * of the configuration, and are not backend specific.
   */
  private def validateExecutableInputs(inputs: ResolvedExecutableInputs, ioFunctions: IoFunctionSet): ErrorOr[Map[OutputPort, WomValue]] = {
    import common.validation.ErrorOr.MapTraversal
    inputs.traverse {
      case (key, value) => value.fold(ResolvedExecutableInputsPoly).apply(ioFunctions) map { key -> _ }
    }
  }

  private def validateWdlFiles(workflowInputs: Map[OutputPort, WomValue]): Checked[Unit] = {

    def prefix(port: OutputPort) = s"Invalid value for File input '${port.fullyQualifiedName}':"

    val failedFiles = workflowInputs collect {
      case (port, womSingleFile: WomSingleFile) if womSingleFile.value.startsWith("\"gs://") =>
        s"""${prefix(port)} ${womSingleFile.value} starts with a '"'"""
      case (port, womSingleFile: WomSingleFile) if womSingleFile.value.isEmpty =>
        s"${prefix(port)} empty value"
    }

    NonEmptyList.fromList(failedFiles.toList) match {
      case Some(errors) => Left(errors)
      case None => Right(())
    }
  }

  def simpleLooksParseable(startsWithOptions: List[String], commentIndicators: List[String])(content: String): Boolean = {
    val fileWithoutInitialWhitespace = content.linesIterator.toList.dropWhile { l =>
      l.forall(_.isWhitespace) || commentIndicators.exists(l.dropWhile(_.isWhitespace).startsWith(_))
    }

    val firstCodeLine = fileWithoutInitialWhitespace.headOption.map(_.dropWhile(_.isWhitespace))
    firstCodeLine.exists { line => startsWithOptions.exists(line.startsWith) }
  }

  def chooseFactory(workflowSource: WorkflowSource, wsfc: WorkflowSourceFilesCollection): ErrorOr[LanguageFactory] = {
    wsfc.workflowType match {
      case Some(languageName) if CromwellLanguages.instance.languages.contains(languageName.toUpperCase) =>
        val language = CromwellLanguages.instance.languages(languageName.toUpperCase)
        wsfc.workflowTypeVersion match {
          case Some(v) if language.allVersions.contains(v) => language.allVersions(v).valid
          case Some(other) => s"Unknown version '$other' for workflow language '$languageName'".invalidNel
          case _ =>
            language.allVersions.values.toList.find(_.looksParsable(workflowSource)).getOrElse(language.default).valid
        }
      case Some(other) => s"Unknown workflow type: $other".invalidNel[LanguageFactory]
      case None =>
        val allFactories = CromwellLanguages.instance.languages.values.flatMap(_.allVersions.values)
        allFactories.find(_.looksParsable(workflowSource)).getOrElse(CromwellLanguages.instance.default.default).validNel
    }
  }

  def findWorkflowSource(workflowSource: Option[WorkflowSource],
                         workflowUrl: Option[WorkflowUrl],
                         resolvers: List[ImportResolver]): Checked[(WorkflowSource, List[ImportResolver])] = {
    (workflowSource, workflowUrl) match {
      case (Some(source), None) => (source, resolvers).validNelCheck
      case (None, Some(url)) =>
        val compoundImportResolver: CheckedAtoB[ImportResolutionRequest, ResolvedImportBundle] = CheckedAtoB.firstSuccess(resolvers.map(_.resolver), s"resolve workflowUrl '$url'")
        val wfSourceAndResolvers: Checked[ResolvedImportBundle] = compoundImportResolver.run(ImportResolutionRequest(url, resolvers))
        wfSourceAndResolvers map { v => (v.source, v.newResolvers) }
      case (Some(_), Some(_)) => "Both workflow source and url can't be supplied".invalidNelCheck
      case (None, None) => "Either workflow source or url has to be supplied".invalidNelCheck
    }
  }
}
