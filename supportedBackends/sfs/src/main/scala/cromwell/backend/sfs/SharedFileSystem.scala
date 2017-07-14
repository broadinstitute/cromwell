package cromwell.backend.sfs

import java.io.{FileNotFoundException, IOException}

import cats.instances.try_._
import cats.syntax.functor._
import com.typesafe.config.Config
import com.typesafe.scalalogging.StrictLogging
import cromwell.backend.io.JobPaths
import cromwell.backend.wdl.WdlFileMapper
import cromwell.core.CromwellFatalExceptionMarker
import cromwell.core.path.{DefaultPath, DefaultPathBuilder, Path, PathFactory}
import lenthall.util.TryUtil
import wdl4s.wdl.EvaluatedTaskInputs
import wdl4s.wdl.values._

import scala.collection.JavaConverters._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object SharedFileSystem extends StrictLogging {

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

  case class PairOfFiles(src: Path, dst: Path)
  type DuplicationStrategy = (Path, Path) => Try[Unit]

  /**
    * Return a `Success` result if the file has already been localized, otherwise `Failure`.
    */
  private def localizePathAlreadyLocalized(originalPath: Path, executionPath: Path): Try[Unit] = {
    if (executionPath.exists) Success(()) else Failure(new RuntimeException(s"$originalPath doesn't exists"))
  }

  private def localizePathViaCopy(originalPath: Path, executionPath: Path): Try[Unit] = {
    val action = Try {
      executionPath.parent.createPermissionedDirectories()
      val executionTmpPath = executionPath.plusExt("tmp")
      originalPath.copyTo(executionTmpPath, overwrite = true).moveTo(executionPath, overwrite = true)
    }.void
    logOnFailure(action, "copy")
  }

  private def localizePathViaHardLink(originalPath: Path, executionPath: Path): Try[Unit] = {
    val action = Try {
      executionPath.parent.createPermissionedDirectories()
      originalPath.linkTo(executionPath)
    }.void
    logOnFailure(action, "hard link")
  }

  private def localizePathViaSymbolicLink(originalPath: Path, executionPath: Path): Try[Unit] = {
      if (originalPath.isDirectory) Failure(new UnsupportedOperationException("Cannot localize directory with symbolic links"))
      else if (!originalPath.exists) Failure(new FileNotFoundException(originalPath.pathAsString))
      else {
        val action = Try {
          executionPath.parent.createPermissionedDirectories()
          executionPath.linkTo(originalPath, symbolic = true)
        }.void
        logOnFailure(action, "symbolic link")
      }
  }

  private def logOnFailure(action: Try[Unit], actionLabel: String): Try[Unit] = {
    if (action.isFailure) logger.warn(s"Localization via $actionLabel has failed: ${action.failed.get.getMessage}")
    action
  }

  private def duplicate(description: String, source: Path, dest: Path, strategies: Stream[DuplicationStrategy]): Try[Unit] = {
    val attempts: Stream[Try[Unit]] = strategies.map(_ (source.followSymbolicLinks, dest))
    attempts.find(_.isSuccess) getOrElse {
      TryUtil.sequence(attempts, s"Could not $description $source -> $dest").void
    }
  }
}

trait SharedFileSystem extends PathFactory {
  import SharedFileSystem._

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

  private def hostAbsoluteFilePath(callRoot: Path, pathString: String): Path = {
    val wdlPath = PathFactory.buildPath(pathString, pathBuilders)
    wdlPath match {
      case _: DefaultPath if !wdlPath.isAbsolute => callRoot.resolve(wdlPath).toAbsolutePath
      case _ => wdlPath
    }
  }

  def outputMapper(job: JobPaths)(wdlValue: WdlValue): Try[WdlValue] = {
    WdlFileMapper.mapWdlFiles(mapJobWdlFile(job))(wdlValue)
  }

  def mapJobWdlFile(job: JobPaths)(wdlFile: WdlFile): WdlFile = {
    wdlFile match {
      case fileNotFound: WdlFile if !hostAbsoluteFilePath(job.callExecutionRoot, fileNotFound.valueString).exists =>
        throw new RuntimeException("Could not process output, file not found: " +
          s"${hostAbsoluteFilePath(job.callExecutionRoot, fileNotFound.valueString).pathAsString}")
      case _ => WdlFile(hostAbsoluteFilePath(job.callExecutionRoot, wdlFile.valueString).pathAsString)
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
    TryUtil.sequenceMap(
      inputs mapValues WdlFileMapper.mapWdlFiles(localizeWdlFile(inputsRoot, docker)),
      "Failures during localization"
    ) recoverWith {
      case e => Failure(new IOException(e.getMessage) with CromwellFatalExceptionMarker)
    }
  }

  def localizeWdlFile(inputsRoot: Path, docker: Boolean)(value: WdlFile): WdlFile = {
    val strategies = if (docker) DockerLocalizers else Localizers

    // Strip the protocol scheme
    def stripProtocolScheme(path: Path): Path = DefaultPathBuilder.get(path.pathWithoutScheme)

    /*
      * Transform an original input path to a path in the call directory.
      * The new path matches the original path, it only "moves" the root to be the call directory.
      */

    def toCallPath(path: String): Try[PairOfFiles] = Try {
      val src = buildPath(path)
      // Strip out potential prefix protocol
      val localInputPath = stripProtocolScheme(src)
      val dest = if (inputsRoot.isParentOf(localInputPath)) localInputPath
      else {
        // Concatenate call directory with absolute input path
        DefaultPathBuilder.get(inputsRoot.pathAsString, localInputPath.pathAsString)
      }

      PairOfFiles(src, dest)
    }

    // Optional function to adjust the path to "docker path" if the call runs in docker
    localizeWdlFile(toCallPath _, strategies.toStream)(value)
  }

  /**
    * Try to localize a WdlFile.
    *
    * @param toDestPath function specifying how to generate the destination path from the source path
    * @param strategies strategies to use for localization
    * @param wdlFile WdlFile to localize
    * @return localized wdl file
    */
  private def localizeWdlFile(toDestPath: (String => Try[PairOfFiles]), strategies: Stream[DuplicationStrategy])
                             (wdlFile: WdlFile): WdlFile = {
    val path = wdlFile.value
    val result = toDestPath(path) flatMap {
      case PairOfFiles(src, dst) => duplicate("localize", src, dst, strategies) map { _ => WdlFile(dst.pathAsString) }
    }
    result.get
  }
}
