package cromwell.languages.util

import cats.data.NonEmptyList
import cats.syntax.validated._
import common.Checked
import common.validation.ErrorOr.ErrorOr
import cromwell.core.CromwellGraphNode.CromwellEnhancedOutputPort
import cromwell.core.path.BetterFileMethods.OpenOptions
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.languages.ValidatedWomNamespace
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
    * @param parentPath if specified, where to unzip to. Otherwise it'll end up somewhere random
    * @return where the imports were unzipped to
    */
  def validateImportsDirectory(zipContents: Array[Byte], parentPath: Option[Path] = None): ErrorOr[Path] = {

    def makeZipFile: Try[Path] = Try {
      DefaultPathBuilder.createTempFile("", ".zip", parentPath).writeByteArray(zipContents)(OpenOptions.default)
    }

    def unZipFile(f: Path) = Try(f.unzipTo(parentPath))

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
}
