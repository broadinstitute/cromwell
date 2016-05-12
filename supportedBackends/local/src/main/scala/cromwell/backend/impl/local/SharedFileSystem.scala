package cromwell.backend.impl.local

import java.nio.file.{FileSystem, Files, Path, Paths}

import com.typesafe.config.Config
import cromwell.backend.BackendJobExecutionActor
import cromwell.core._
import wdl4s.CallInputs
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.types.{WdlArrayType, WdlMapType}
import wdl4s.util.TryUtil
import wdl4s.values.{WdlValue, _}

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

trait SharedFileSystem extends PathFactory { this: BackendJobExecutionActor =>
  import SharedFileSystem._
  import better.files._

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

  def processOutputs(engineFunctions: WdlStandardLibraryFunctions, jobPaths: JobPaths): Try[CallOutputs] = {
    val evaluatedOutputs = evaluateOutputs(engineFunctions, outputMapper(jobPaths))
    lazy val workflowId = jobDescriptor.descriptor.id

    TryUtil.sequenceMap(evaluatedOutputs, s"Workflow $workflowId post processing failed")
  }

  private def hostAbsoluteFilePath(callRoot: Path, pathString: String): Path = {
    val wdlPath = Paths.get(pathString)
    callRoot.resolve(wdlPath).toAbsolutePath
  }

  private def outputMapper(job: JobPaths)(wdlValue: WdlValue): Try[WdlValue] = {
    wdlValue match {
      case fileNotFound: WdlFile if !hostAbsoluteFilePath(job.callRoot, fileNotFound.valueString).exists =>
        Failure(new RuntimeException(s"Could not process output, file not found: ${hostAbsoluteFilePath(job.callRoot, fileNotFound.valueString).toString}"))
      case file: WdlFile => Try(WdlFile(hostAbsoluteFilePath(job.callRoot, file.valueString).toString))
      case array: WdlArray =>
        val mappedArray = array.value map outputMapper(job)
        TryUtil.sequence(mappedArray) map { WdlArray(array.wdlType, _) }
      case other => Success(other)
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
      val src = buildPath(path, filesystems)
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
}
