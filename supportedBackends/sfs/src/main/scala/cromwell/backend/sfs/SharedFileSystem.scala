package cromwell.backend.sfs

import java.io.{FileNotFoundException, IOException}
import java.nio.file.{Path, Paths}

import cats.instances.try_._
import cats.syntax.functor._
import com.typesafe.config.Config
import com.typesafe.scalalogging.StrictLogging
import cromwell.backend.io.JobPaths
import cromwell.core.CromwellFatalExceptionMarker
import cromwell.core.path.PathFactory
import cromwell.core.path.PathFactory._
import lenthall.util.TryUtil
import wdl4s.EvaluatedTaskInputs
import wdl4s.types.{WdlArrayType, WdlMapType}
import wdl4s.values._

import scala.collection.JavaConverters._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object SharedFileSystem extends StrictLogging {
  import better.files._

  final case class AttemptedLookupResult(name: String, value: Try[WdlValue]) {
    def toPair: (String, Try[WdlValue]) = name -> value
  }

  object AttemptedLookupResult {
    implicit class AugmentedAttemptedLookupSequence(s: Seq[AttemptedLookupResult]) {
      def toLookupMap: Map[String, WdlValue] = s collect {
        case AttemptedLookupResult(name, Success(value)) => (name, value)
      } toMap
    }
  }

  case class PairOfFiles(src: File, dst: File)
  type DuplicationStrategy = (File, File) => Try[Unit]

  /**
    * Return a `Success` result if the file has already been localized, otherwise `Failure`.
    */
  private def localizePathAlreadyLocalized(originalPath: File, executionPath: File): Try[Unit] = {
    if (executionPath.exists) Success(()) else Failure(new RuntimeException(s"$originalPath doesn't exists"))
  }

  private def localizePathViaCopy(originalPath: File, executionPath: File): Try[Unit] = {
    val action = Try {
      executionPath.parent.createDirectories()
      val executionTmpPath = pathPlusSuffix(executionPath, "tmp")
      originalPath.copyTo(executionTmpPath, overwrite = true).moveTo(executionPath, overwrite = true)
    }.void
    logOnFailure(action, "copy")
  }

  private def localizePathViaHardLink(originalPath: File, executionPath: File): Try[Unit] = {
    val action = Try {
      executionPath.parent.createDirectories()
      executionPath.linkTo(originalPath, symbolic = false)
    }.void
    logOnFailure(action, "hard link")
  }

  private def localizePathViaSymbolicLink(originalPath: File, executionPath: File): Try[Unit] = {
      if (originalPath.isDirectory) Failure(new UnsupportedOperationException("Cannot localize directory with symbolic links"))
      else if (!originalPath.exists) Failure(new FileNotFoundException(originalPath.pathAsString))
      else {
        val action = Try {
          executionPath.parent.createDirectories()
          executionPath.linkTo(originalPath, symbolic = true)
        }.void
        logOnFailure(action, "symbolic link")
      }
  }

  private def logOnFailure(action: Try[Unit], actionLabel: String): Try[Unit] = {
    if (action.isFailure) logger.warn(s"Localization via $actionLabel has failed: ${action.failed.get.getMessage}")
    action
  }

  private def duplicate(description: String, source: File, dest: File, strategies: Stream[DuplicationStrategy]): Try[Unit] = {
    import cromwell.util.FileUtil._

    val attempts: Stream[Try[Unit]] = strategies.map(_ (source.followSymlinks, dest))
    attempts.find(_.isSuccess) getOrElse {
      TryUtil.sequence(attempts, s"Could not $description $source -> $dest").void
    }
  }
}

trait SharedFileSystem extends PathFactory {
  import SharedFileSystem._
  import better.files._

  def sharedFileSystemConfig: Config

  lazy val DefaultStrategies = Seq("hard-link", "soft-link", "copy")

  lazy val LocalizationStrategies: Seq[String] = getConfigStrategies("localization")
  lazy val Localizers: Seq[DuplicationStrategy] = createStrategies(LocalizationStrategies, docker = false)
  lazy val DockerLocalizers: Seq[DuplicationStrategy] = createStrategies(LocalizationStrategies, docker = true)

  lazy val CachingStrategies: Seq[String] = getConfigStrategies("caching.duplication-strategy")
  lazy val Cachers: Seq[DuplicationStrategy] = createStrategies(CachingStrategies, docker = false)

  private def getConfigStrategies(configPath: String): Seq[String] = {
    if (sharedFileSystemConfig.hasPath(configPath)) {
      sharedFileSystemConfig.getStringList(configPath).asScala
    } else {
      DefaultStrategies
    }
  }

  private def createStrategies(configStrategies: Seq[String], docker: Boolean): Seq[DuplicationStrategy] = {
    // If localizing for a docker job, remove soft-link as an option
    val filteredConfigStrategies = configStrategies filter {
      case "soft-link" if docker => false
      case _ => true
    }

    // Convert the (remaining) config strategies to duplication strategies
    val mappedDuplicationStrategies = filteredConfigStrategies map {
      case "hard-link" => localizePathViaHardLink _
      case "soft-link" => localizePathViaSymbolicLink _
      case "copy" => localizePathViaCopy _
      case unsupported => throw new UnsupportedOperationException(s"Strategy $unsupported is not recognized")
    }

    // Prepend the default duplication strategy, and return the sequence
    localizePathAlreadyLocalized _ +: mappedDuplicationStrategies
  }

  private def hostAbsoluteFilePath(callRoot: Path, pathString: String): File = {
    val wdlPath = PathFactory.buildPath(pathString, pathBuilders)
    callRoot.resolve(wdlPath).toAbsolutePath
  }

  def outputMapper(job: JobPaths)(wdlValue: WdlValue): Try[WdlValue] = {
    wdlValue match {
      case fileNotFound: WdlFile if !hostAbsoluteFilePath(job.callExecutionRoot, fileNotFound.valueString).exists =>
        Failure(new RuntimeException("Could not process output, file not found: " +
          s"${hostAbsoluteFilePath(job.callExecutionRoot, fileNotFound.valueString).pathAsString}"))
      case file: WdlFile => Try(WdlFile(hostAbsoluteFilePath(job.callExecutionRoot, file.valueString).pathAsString))
      case array: WdlArray =>
        val mappedArray = array.value map outputMapper(job)
        TryUtil.sequence(mappedArray) map { WdlArray(array.wdlType, _) }
      case map: WdlMap =>
        val mappedMap = map.value mapValues outputMapper(job)
        TryUtil.sequenceMap(mappedMap) map { WdlMap(map.wdlType, _) }
      case other => Success(other)
    }
  }

  def cacheCopy(sourceFilePath: Path, destinationFilePath: Path): Try[Unit] = {
    duplicate("cache", sourceFilePath, destinationFilePath, Cachers.toStream)
  }

  /**
   * Return a possibly altered copy of inputs reflecting any localization of input file paths that might have
   * been performed for this `Backend` implementation.
   */
  def localizeInputs(inputsRoot: Path, docker: Boolean)(inputs: EvaluatedTaskInputs): Try[EvaluatedTaskInputs] = {
    val strategies = if (docker) DockerLocalizers else Localizers

    // Use URI to identify protocol scheme and strip it out
    def stripProtocolScheme(path: Path): Path = {
      val uri = path.toUri
      val host = Option(uri.getHost)
      val uriPath = uri.getPath

      host map { h => Paths.get(h, uriPath) } getOrElse Paths.get(uriPath)
    }

    /*
      * Transform an original input path to a path in the call directory.
      * The new path matches the original path, it only "moves" the root to be the call directory.
      */

    def toCallPath(path: String): Try[PairOfFiles] = Try {
      val src = buildFile(path)
      // Strip out potential prefix protocol
      val localInputPath = stripProtocolScheme(src.path)
      val dest = if (File(inputsRoot).isParentOf(localInputPath)) File(localInputPath)
      else {
        // Concatenate call directory with absolute input path
        File(Paths.get(inputsRoot.toString, localInputPath.toString))
      }

      PairOfFiles(src, dest)
    }

    // Optional function to adjust the path to "docker path" if the call runs in docker
    val localizeFunction = localizeWdlValue(toCallPath, strategies.toStream) _
    val localizedValues = inputs.toSeq map {
      case (declaration, value) => localizeFunction(value) map { declaration -> _ }
    }

    TryUtil.sequence(localizedValues, "Failures during localization").map(_.toMap) recoverWith {
      case e => Failure(new IOException(e.getMessage) with CromwellFatalExceptionMarker)
    }
  }

  /**
   * Try to localize a WdlValue if it is or contains a WdlFile.
   *
   * @param toDestPath function specifying how to generate the destination path from the source path
   * @param strategies strategies to use for localization
   * @param wdlValue WdlValue to localize
   * @return localized wdlValue
   */
  private def localizeWdlValue(toDestPath: (String => Try[PairOfFiles]), strategies: Stream[DuplicationStrategy])
                              (wdlValue: WdlValue): Try[WdlValue] = {

    def adjustArray(t: WdlArrayType, inputArray: Seq[WdlValue]): Try[WdlArray] = {
      val tryAdjust = inputArray map localizeWdlValue(toDestPath, strategies)

      TryUtil.sequence(tryAdjust, s"Failed to localize files in input Array ${wdlValue.valueString}") map { adjusted =>
        new WdlArray(t, adjusted)
      }
    }

    def adjustMap(t: WdlMapType, inputMap: Map[WdlValue, WdlValue]): Try[WdlMap] = {
      val tryAdjust = inputMap mapValues { localizeWdlValue(toDestPath, strategies) }

      TryUtil.sequenceMap(tryAdjust, s"Failed to localize files in input Map ${wdlValue.valueString}") map { adjusted =>
        new WdlMap(t, adjusted)
      }
    }

    def adjustFile(path: String) = {
      toDestPath(path) flatMap {
        case PairOfFiles(src, dst) => duplicate("localize", src, dst, strategies) map { _ => WdlFile(dst.toString) }
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
