package cromwell.backend.impl.local

import java.io.File
import java.nio.file.{FileSystem, Files, Path, Paths}

import com.typesafe.config.Config
import cromwell.backend.BackendJobDescriptor
import cromwell.core._
import wdl4s.WdlExpression.ScopedLookupFunction
import wdl4s.expression.WdlFunctions
import wdl4s.types.{WdlArrayType, WdlFileType, WdlMapType}
import wdl4s.util.TryUtil
import wdl4s.values.{WdlValue, _}
import wdl4s.{CallInputs, TaskOutput}

import scala.collection.JavaConverters._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object SharedFileSystem {
  import better.files._

  final case class AttemptedLookupResult(name: String, value: Try[WdlValue]) {
    def toPair = name -> value
  }

  object AttemptedLookupResult {
    implicit class AugmentedAttemptedLookupSequence(s: Seq[AttemptedLookupResult]) {
      def toLookupMap: Map[String, WdlValue] = s collect {
        case AttemptedLookupResult(name, Success(value)) => (name, value)
      } toMap
    }
  }

  type PathsPair = (Path, Path)
  type LocalizationStrategy = (Path, Path) => Try[Unit]

  /**
    * Return a `Success` result if the file has already been localized, otherwise `Failure`.
    */
  private def localizePathAlreadyLocalized(originalPath: Path, executionPath: Path): Try[Unit] = {
    if (executionPath.exists) Success(()) else Failure(new Throwable)
  }

  private def localizePathViaCopy(originalPath: Path, executionPath: Path): Try[Unit] = {
    executionPath.getParent.createDirectories()
    Try(originalPath.copyTo(executionPath))
  }

  private def localizePathViaHardLink(originalPath: Path, executionPath: Path): Try[Unit] = {
    executionPath.getParent.createDirectories()
    Try(Files.createLink(executionPath, originalPath))
  }

  /**
    * TODO: The 'call' parameter here represents the call statement in WDL that references this path.
    * We're supposed to not use symbolic links if the call uses Docker.  However, this is currently a
    * bit incorrect because multiple calls can reference the same path if that path is in a declaration.
    *
    * The symbolic link will only fail in the Docker case if a Call uses the file directly and not
    * indirectly through one of its input expressions
    */

  private def localizePathViaSymbolicLink(originalPath: Path, executionPath: Path): Try[Unit] = {
      if (originalPath.isDirectory) Failure(new UnsupportedOperationException("Cannot localize directory with symbolic links"))
      else {
        executionPath.getParent.createDirectories()
        Try(Files.createSymbolicLink(executionPath, originalPath.toAbsolutePath))
      }
  }
}

trait SharedFileSystem {
  import SharedFileSystem._

  def sharedFsConfig: Config

  lazy val LocalizationStrategies = sharedFsConfig.getStringList("localization").asScala
  lazy val Localizers = localizePathAlreadyLocalized _ +: (LocalizationStrategies map {
    case "hard-link" => localizePathViaHardLink _
    case "soft-link" => localizePathViaSymbolicLink _
    case "copy" => localizePathViaCopy _
    case unsupported => throw new UnsupportedOperationException(s"Localization strategy $unsupported is not recognized")
  })

  // Note that any unrecognized configuration will be raised when Localizers (just above) gets resolved.
  lazy val DockerLocalizers = localizePathAlreadyLocalized _ +: (LocalizationStrategies collect {
    case "hard-link" => localizePathViaHardLink _
    case "copy" => localizePathViaCopy _
  })

  def processOutputs(jobDescriptor: BackendJobDescriptor, workflowId: WorkflowId, lookup: ScopedLookupFunction, engineFunctions: WdlFunctions[WdlValue], jobPaths: JobPaths): Try[CallOutputs] = {
    def outputFoldingFunction = {
      (currentList: Seq[AttemptedLookupResult], taskOutput: TaskOutput) => {
        currentList ++ Seq(AttemptedLookupResult(taskOutput.name, outputLookup(taskOutput, currentList)))
      }
    }

    def outputLookup(taskOutput: TaskOutput, currentList: Seq[AttemptedLookupResult]) = for {
      expressionValue <- taskOutput.requiredExpression.evaluate(lookup, engineFunctions)
      convertedValue <- outputAutoConversion(jobDescriptor, taskOutput, expressionValue, jobPaths)
      pathAdjustedValue <- Success(absolutizeOutputWdlFile(convertedValue, jobPaths.callRoot))
    } yield pathAdjustedValue

    val outputs = jobDescriptor.key.call.task.outputs
    val outputMappings = outputs.foldLeft(Seq.empty[AttemptedLookupResult])(outputFoldingFunction).map(_.toPair).toMap

    TryUtil.sequenceMap(outputMappings, s"Workflow $workflowId post processing failed") map {
      _.mapValues { wdlValue => CallOutput(wdlValue, None) }
    }
  }
  /**
   * Return a possibly altered copy of inputs reflecting any localization of input file paths that might have
   * been performed for this `Backend` implementation.
   * NOTE: This ends up being a backdoor implementation of Backend.adjustInputPaths as both LocalBackend and SgeBackend
   *    end up with this implementation and thus use it to satisfy their contract with Backend.
   *    This is yuck-tastic and I consider this a FIXME, but not for this refactor
   */
  def localizeInputs(jobPaths: JobPaths, docker: Boolean, filesystems: List[FileSystem], inputs: CallInputs): CallInputs = {

    val strategies = if (docker) DockerLocalizers else Localizers

    // Use URI to identify protocol scheme and strip it out
    def stripProtocolScheme(path: Path): Path = {
      val uri = path.toUri
      val host = Option(uri.getHost)
      val uriPath = uri.getPath

      host map { h => Paths.get(h).resolve(uriPath) } getOrElse Paths.get(uriPath)
    }

    /**
      * Transform an original input path to a path in the call directory.
      * The new path matches the original path, it only "moves" the root to be the call directory.
      */
    def toCallPath(path: String): PathsPair = {
      val src = PathFactory.buildPath(path, filesystems)
      // Strip out potential prefix protocol
      val localInputPath = stripProtocolScheme(src).toString
      // Concatenate call directory with absolute input path
      (src, Paths.get(jobPaths.callRoot.toAbsolutePath.toString, localInputPath.toString))
    }

    // Optional function to adjust the path to "docker path" if the call runs in docker
    val postProcessor: Option[Path => Path] = if (docker) Option(jobPaths.toDockerPath _) else None
    val localizeFunction = localizeWdlValue(toCallPath, strategies.toStream, postProcessor) _
    val localizedValues = inputs.toSeq map {
      case (name, value) => localizeFunction(value) map { name -> _ }
    }

    TryUtil.sequence(localizedValues, "Failures during localization").get toMap
  }

  /**
   * Try to localize a WdlValue if it is or contains a WdlFile.
   *
   * @param toDestPath function specifying how to generate the destination path from the source path
   * @param strategies strategies to use for localization
   * @param postProcessor optional function to be applied to the path after the file it points to has been localized (defaults to noop)
   * @param wdlValue WdlValue to localize
   * @return localized wdlValue
   */
  def localizeWdlValue(toDestPath: (String => PathsPair), strategies: Stream[LocalizationStrategy], postProcessor: Option[Path => Path] = None)(wdlValue: WdlValue): Try[WdlValue] = {

    def localize(source: Path, dest: Path) = strategies map { _(source, dest) } find { _.isSuccess } getOrElse {
      Failure(new UnsupportedOperationException(s"Could not localize $source -> $dest"))
    }

    def adjustArray(t: WdlArrayType, inputArray: Seq[WdlValue]): Try[WdlArray] = {
      val tryAdjust = inputArray map localizeWdlValue(toDestPath, strategies, postProcessor)

      TryUtil.sequence(tryAdjust, s"Failed to localize files in input Array ${wdlValue.valueString}") map { adjusted =>
        new WdlArray(t, adjusted)
      }
    }

    def adjustMap(t: WdlMapType, inputMap: Map[WdlValue, WdlValue]): Try[WdlMap] = {
      val tryAdjust = inputMap mapValues { localizeWdlValue(toDestPath, strategies, postProcessor) }

      TryUtil.sequenceMap(tryAdjust, s"Failed to localize files in input Map ${wdlValue.valueString}") map { adjusted =>
        new WdlMap(t, adjusted)
      }
    }

    def adjustFile(path: String) = {
      val (src, dst) = toDestPath(path)
      localize(src, dst) map { Unit =>
        val finalPath = postProcessor map { _(dst) } getOrElse dst
        WdlFile(finalPath.toAbsolutePath.toString)
      }
    }

    wdlValue match {
      case wdlFile: WdlFile => adjustFile(wdlFile.value)
      case WdlArray(t, values) => adjustArray(t, values)
      case WdlMap(t, values) => adjustMap(t, values)
      case x => Success(x)
    }
  }

  private def assertTaskOutputPathExists(path: String, taskOutput: TaskOutput, callFqn: String): Try[WdlFile] =
    if (Files.exists(Paths.get(path))) Success(WdlFile(path))
    else Failure(new RuntimeException(
      s"""ERROR: Could not process output '${taskOutput.wdlType.toWdlString} ${taskOutput.name}' of $callFqn:
         |
         |Invalid path: $path
       """.stripMargin
    ))

  private def hostAbsoluteFilePath(callRoot: Path, pathString: String): String =
    if (new File(pathString).isAbsolute) pathString else Paths.get(callRoot.toAbsolutePath.toString, pathString).toString

  /**
   * Handle possible auto-conversion from an output expression `WdlString` to a `WdlFile` task output.
   * The following should work:
   *
   * <pre>
   * File bam = "foo.bam"
   * </pre>
   *
   * This also handles coercions to target types.  For example, a task output may look like this:
   *
   * <pre>
   * Map[Int, String] my_map = read_map(stdout())
   * </pre>
   *
   * read_map() will return a Map[String, String] but if the target type is Map[Int, String], this
   * function will attempt to do the coercion.
   *
   * read_lines() will return an Array[String] but if the target type is Array[Float], then this
   * function will do that conversion.
   */
  private def outputAutoConversion(jobDescriptor: BackendJobDescriptor, taskOutput: TaskOutput, rawOutputValue: WdlValue, job: JobPaths): Try[WdlValue] = {
    rawOutputValue match {
      case rhs if rhs.wdlType == taskOutput.wdlType => Success(rhs)
      case rhs: WdlString if taskOutput.wdlType == WdlFileType => assertTaskOutputPathExists(hostAbsoluteFilePath(job.callRoot, rhs.value), taskOutput, jobDescriptor.key.call.fullyQualifiedName)
      case rhs => taskOutput.wdlType.coerceRawValue(rhs)
    }
  }

  private def absolutizeOutputWdlFile(value: WdlValue, cwd: Path): WdlValue = value match {
    case f: WdlFile if f.valueString(0) != '/' => WdlFile(cwd.resolve(f.valueString).toAbsolutePath.toString)
    case a: WdlArray => a map {absolutizeOutputWdlFile(_, cwd)}
    case x => x
  }
}
