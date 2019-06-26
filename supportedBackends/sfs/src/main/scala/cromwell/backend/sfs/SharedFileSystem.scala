package cromwell.backend.sfs

import java.io.{FileNotFoundException, IOException}

import akka.actor.ActorContext
import akka.stream.ActorMaterializer
import cats.instances.try_._
import cats.syntax.functor._
import com.typesafe.config.Config
import com.typesafe.scalalogging.StrictLogging
import common.collections.EnhancedCollections._
import common.util.TryUtil
import cromwell.backend.io.JobPaths
import cromwell.core.CromwellFatalExceptionMarker
import cromwell.core.path.{DefaultPath, DefaultPathBuilder, Path, PathFactory}
import cromwell.filesystems.http.HttpPathBuilder
import net.ceedubs.ficus.Ficus._
import wom.WomFileMapper
import wom.values._

import scala.collection.JavaConverters._
import scala.collection.mutable
import scala.concurrent.Await
import scala.concurrent.duration.Duration
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object SharedFileSystem extends StrictLogging {

  final case class AttemptedLookupResult(name: String, value: Try[WomValue]) {
    def toPair: (String, Try[WomValue]) = name -> value
  }

  object AttemptedLookupResult {
    implicit class AugmentedAttemptedLookupSequence(s: Seq[AttemptedLookupResult]) {
      def toLookupMap: Map[String, WomValue] = s collect {
        case AttemptedLookupResult(name, Success(value)) => (name, value)
      } toMap
    }
  }

  case class PairOfFiles(src: Path, dst: Path)
  type DuplicationStrategy = (Path, Path, Boolean) => Try[Unit]
  
  private def createParentDirectory(executionPath: Path, docker: Boolean) = {
    if (docker) executionPath.parent.createPermissionedDirectories()
    else executionPath.parent.createDirectories()
  }

  /**
    * Return a `Success` result if the file has already been localized, otherwise `Failure`.
    */
  private def localizePathAlreadyLocalized(originalPath: Path, executionPath: Path, docker: Boolean): Try[Unit] = {
    if (executionPath.exists) Success(()) else Failure(new RuntimeException(s"$originalPath doesn't exist"))
  }

  private def localizePathViaCopy(originalPath: Path, executionPath: Path, docker: Boolean): Try[Unit] = {
    val action = Try {
      createParentDirectory(executionPath, docker)
      val executionTmpPath = executionPath.plusExt("tmp")
      originalPath.copyTo(executionTmpPath, overwrite = true).moveTo(executionPath, overwrite = true)
    }.void
    logOnFailure(action, "copy")
  }

  private def localizePathViaHardLink(originalPath: Path, executionPath: Path, docker: Boolean): Try[Unit] = {
    val action = Try {
      createParentDirectory(executionPath, docker)
      originalPath.linkTo(executionPath)
    }.void
    logOnFailure(action, "hard link")
  }

  private def localizePathViaSymbolicLink(originalPath: Path, executionPath: Path, docker: Boolean): Try[Unit] = {
      if (originalPath.isDirectory) Failure(new UnsupportedOperationException("Cannot localize directory with symbolic links"))
      else if (!originalPath.exists) Failure(new FileNotFoundException(originalPath.pathAsString))
      else {
        val action = Try {
          createParentDirectory(executionPath, docker)
          executionPath.linkTo(originalPath, symbolic = true)
        }.void
        logOnFailure(action, "symbolic link")
      }
  }

  private def logOnFailure(action: Try[Unit], actionLabel: String): Try[Unit] = {
    if (action.isFailure) logger.warn(s"Localization via $actionLabel has failed: ${action.failed.get.getMessage}")
    action
  }

  private def duplicate(description: String, source: Path, dest: Path, strategies: Stream[DuplicationStrategy], docker: Boolean): Try[Unit] = {
    val attempts: Stream[Try[Unit]] = strategies.map(_.apply(source.followSymbolicLinks, dest, docker))
    attempts.find(_.isSuccess) getOrElse {
      TryUtil.sequence(attempts, s"Could not $description $source -> $dest").void
    }
  }

  private lazy val beingCopied: mutable.Map[Path, Boolean] = mutable.Map[Path, Boolean]()

  private def waitOnCopy(path: Path, lockFile: Path): Unit = {
    while (beingCopied.getOrElse(path, false) || lockFile.exists)  {
      Thread.sleep(1)
    }
  }

  private def countLinks(path: Path): Int = {
    path.getAttribute("unix:nlink").asInstanceOf[Int]
  }
}

trait SharedFileSystem extends PathFactory {
  import SharedFileSystem._

  def sharedFileSystemConfig: Config
  lazy val maxHardLinks: Int = sharedFileSystemConfig.getOrElse[Int]("max-hardlinks",950)  // Windows limit 1024. Keep a safe margin.
  lazy val cachedCopyDir: Option[Path] = None

  private def localizePathViaCachedCopy(originalPath: Path, executionPath: Path, docker: Boolean): Try[Unit] = {
    val action = Try {
      createParentDirectory(executionPath, docker)
      // Hash the parent. This will make sure bamfiles and their indexes stay in the same dir. This is not ideal. But should work.
      // There should be no file collisions because two files with the same name cannot exist in the same parent dir.
      // use .get . This strategy should not be used when there is no cachedCopyDir
      val cachedCopySubDir: Path = cachedCopyDir.get.createChild(originalPath.toAbsolutePath.parent.hashCode.toString, asDirectory = true)

      // By prepending the modtime we prevent collisions in the cache from files that have changed in between.
      // Md5 is safer but much much slower and way too CPU intensive for big files.
      val pathAndModTime: String = originalPath.lastModifiedTime.toEpochMilli.toString + originalPath.name
      val cachedCopyPath: Path = cachedCopySubDir./(pathAndModTime)
      val cachedCopyPathLockFile: Path = cachedCopyPath.plusSuffix(".lock")

      if (!cachedCopyPath.exists || countLinks(cachedCopyPath) >= maxHardLinks) {
        // This variable is used so we can release the lock before we start with the copying.
        var shouldCopy = false

        // LOCK: This block should only be accessed by one thread at a time. Otherwise multiple threads
        // will copy the same cache file. The performance impact should be minimal as it is only a few
        // decisions which are not time consuming.
        SharedFileSystem.synchronized {

          // We check again if cachedCopyPath is there or if the number of Hardlinks is still exceeded.
          // The copying may have been started while waiting on the lock.
          // If it is not there or the maxHardLinks are exceeded, is it already being copied by another thread?
          // if not copied by another thread, is it copied by another cromwell process? (Lock file present)
          if ((!cachedCopyPath.exists || countLinks(cachedCopyPath) >= maxHardLinks) &&
            !SharedFileSystem.beingCopied.getOrElse(cachedCopyPath, false) &&
            !cachedCopyPathLockFile.exists) {
            // Create a lock file so other cromwell processes know copying has started
            try {
              cachedCopyPathLockFile.touch()
              // Create an entry in the dictionary so other threads in this process know copying
              // has started (filesystem can be too slow if multiple threads need the same file at
              // exactly the same time)
              SharedFileSystem.beingCopied(cachedCopyPath) = true
              // Set should copy to true for *this* thread. The lock can now be released.
              // this thread will do the copying. The other threads and/or processes will wait for the copying to
              // be completed.
              shouldCopy = true
            } catch {
              case e: Exception =>
                // In case any error happens, make sure that all locks are removed. Otherwise cromwell will hang in the
                // future.
                SharedFileSystem.beingCopied.remove(cachedCopyPath)
                cachedCopyPathLockFile.delete(true)
                throw e
            }
          }
        } // LOCK RELEASE
        if (shouldCopy) {
          try {
            val cachedCopyTmpPath = cachedCopyPath.plusExt("tmp")
            // CachedCopyPath is overwritten. It is possible that the number of hardlinks is exceeded. In which case
            // the file is already there.
            originalPath.copyTo(cachedCopyTmpPath, overwrite = true).moveTo(cachedCopyPath, overwrite = true)
          } catch {
            case e: Exception => throw e
          }
          finally {
            // Always remove the locks after copying. Even if there is an exception.
            // We remove the key! Not set it to false. We don't want this map being flooded with
            // keys if the cromwell process is active for months in server mode. (Memory leak!)
            SharedFileSystem.beingCopied.remove(cachedCopyPath)
            cachedCopyPathLockFile.delete()
          }
        }
      }
      // If the file does exist, it may still be copied/moved and therefore not all the bits might
      // be there. Therefore we always wait until the locks are cleared.
      SharedFileSystem.waitOnCopy(cachedCopyPath, cachedCopyPathLockFile)
      cachedCopyPath.linkTo(executionPath)
    }.void
    logOnFailure(action, "cached copy file")
  }

  implicit def actorContext: ActorContext

  lazy val DefaultStrategies = Seq("hard-link", "soft-link", "copy")

  lazy val LocalizationStrategyNames: Seq[String] = getConfigStrategies("localization")
  lazy val LocalizationStrategies: Seq[DuplicationStrategy] = createStrategies(LocalizationStrategyNames, docker = false)
  lazy val DockerLocalizationStrategies: Seq[DuplicationStrategy] = createStrategies(LocalizationStrategyNames, docker = true)

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
    // If no cachedCopyDir is defined, cached-copy can not be used and is removed.
    val filteredConfigStrategies = configStrategies filter {
      case "soft-link" if docker => false
      case "cached-copy" if cachedCopyDir.isEmpty => false
      case _ => true
    }

    // Convert the (remaining) config strategies to duplication strategies
    val mappedDuplicationStrategies = filteredConfigStrategies map {
      case "hard-link" => localizePathViaHardLink _
      case "soft-link" => localizePathViaSymbolicLink _
      case "copy" => localizePathViaCopy _
      case "cached-copy" => localizePathViaCachedCopy _
      case unsupported => throw new UnsupportedOperationException(s"Strategy $unsupported is not recognized")
    }

    // Prepend the default duplication strategy, and return the sequence
    localizePathAlreadyLocalized _ +: mappedDuplicationStrategies
  }

  def hostAbsoluteFilePath(jobPaths: JobPaths, pathString: String): Path = {
    val path = PathFactory.buildPath(pathString, pathBuilders)
    path match {
      case _: DefaultPath if !path.isAbsolute => jobPaths.callExecutionRoot.resolve(path).toAbsolutePath
      case _: DefaultPath if jobPaths.isInExecution(path.pathAsString) => jobPaths.hostPathFromContainerPath(path.pathAsString)
      case _: DefaultPath => jobPaths.hostPathFromContainerInputs(path.pathAsString)
    }
  }

  def outputMapper(job: JobPaths)(womValue: WomValue): Try[WomValue] = {
    WomFileMapper.mapWomFiles(mapJobWomFile(job), Set.empty)(womValue)
  }

  def mapJobWomFile(jobPaths: JobPaths)(womFile: WomFile): WomFile = {
    val hostPath = hostAbsoluteFilePath(jobPaths, womFile.valueString)
    def hostAbsolute(pathString: String): String = hostAbsoluteFilePath(jobPaths, pathString).pathAsString

    if (!hostPath.exists) throw new FileNotFoundException(s"Could not process output, file not found: ${hostAbsolute(womFile.valueString)}")

    // There are composite WomFile types like WomMaybeListedDirectoryType that need to make the paths of contained
    // WomFiles host absolute, so don't just pass in a `const` of the function call result above.
    womFile mapFile hostAbsolute
  }

  def cacheCopy(sourceFilePath: Path, destinationFilePath: Path): Try[Unit] = {
    duplicate("cache", sourceFilePath, destinationFilePath, Cachers.toStream, docker = false)
  }

  /**
   * Return a possibly altered copy of inputs reflecting any localization of input file paths that might have
   * been performed for this `Backend` implementation.
   */
  def localizeInputs(inputsRoot: Path, docker: Boolean)(inputs: WomEvaluatedCallInputs): Try[WomEvaluatedCallInputs] = {
    TryUtil.sequenceMap(
      inputs safeMapValues WomFileMapper.mapWomFiles(localizeWomFile(inputsRoot, docker), Set.empty),
      "Failures during localization"
    ) recoverWith {
      case e => Failure(new IOException(e.getMessage) with CromwellFatalExceptionMarker)
    }
  }

  def localizeWomFile(inputsRoot: Path, docker: Boolean)(value: WomFile): WomFile = {
    val strategies = if (docker) DockerLocalizationStrategies else LocalizationStrategies

    // Strip the protocol scheme
    def stripProtocolScheme(path: Path): Path = DefaultPathBuilder.get(path.pathWithoutScheme)

    /*
      * Transform an original input path to a path in the call directory.
      * The new path matches the original path, it only "moves" the root to be the call directory.
      */

    def toCallPath(womFile: WomFile)(path: String): Try[PairOfFiles] = Try {
      val src = buildPath(path)
      // Strip out potential prefix protocol
      val localInputPath = stripProtocolScheme(src)
      val dest = if (inputsRoot.isParentOf(localInputPath)) localInputPath
      else {
        val nameOverride = womFile match {
          case directory: WomMaybeListedDirectory => directory.basename
          case _ => None
        }
        // Concatenate call directory with absolute input path
        DefaultPathBuilder.get(inputsRoot.pathAsString, localInputPath.parent.pathAsString.hashCode.toString, nameOverride.getOrElse(localInputPath.name))
      }

      PairOfFiles(src, dest)
    }

    // A possibly staged version of the input file suitable for downstream processing, or just the original input
    // file if no staging was required.
    val staged: WomFile = value.mapFile { input =>
      pathBuilders.collectFirst({ case h: HttpPathBuilder if HttpPathBuilder.accepts(input) => h }) match {
        case Some(httpPathBuilder) =>
          implicit val materializer = ActorMaterializer()
          implicit val ec = actorContext.dispatcher

          Await.result(httpPathBuilder.content(input).map { _.toString }, Duration.Inf)
        case _ => input
      }
    }

    // Optional function to adjust the path to "docker path" if the call runs in docker
    localizeWomFile(toCallPath _, strategies.toStream, docker)(staged)
  }

  /**
    * Try to localize a WomFile.
    *
    * @param toDestPath function specifying how to generate the destination path from the source path
    * @param strategies strategies to use for localization
    * @param womFile WomFile to localize
    * @return localized WomFile
    */
  private def localizeWomFile(toDestPath: WomFile => String => Try[PairOfFiles], strategies: Stream[DuplicationStrategy], docker: Boolean)
                             (womFile: WomFile): WomFile = {
    val localized = womFile mapWomFile { file =>
      val result = toDestPath(file)(file.value) flatMap {
        case PairOfFiles(src, dst) => duplicate("localize", src, dst, strategies, docker).map(_ => dst.pathAsString)
      }
      result.get
    }
    val sized = localized collect {
      case womMaybePopulatedFile@WomMaybePopulatedFile(Some(path), _, None, _, _, _, _) =>
        val pair = toDestPath(womMaybePopulatedFile)(path).get
        val srcSize = pair.src.size
        womMaybePopulatedFile.copy(sizeOption = Option(srcSize))
    }
    sized
  }
}
