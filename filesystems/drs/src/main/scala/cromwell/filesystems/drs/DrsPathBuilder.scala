package cromwell.filesystems.drs

import cloud.nio.impl.drs.DrsCloudNioFileSystemProvider
import com.typesafe.scalalogging.StrictLogging
import cromwell.core.path.PathFactory.PathBuilders
import cromwell.core.path.{Path, PathFactory, PreResolvePathBuilder}

import scala.util.{Failure, Success, Try}

case class DrsPathBuilder(fileSystemProvider: DrsCloudNioFileSystemProvider,
                          requesterPaysProjectIdOption: Option[String],
                          preResolve: Boolean = false,
                         ) extends PreResolvePathBuilder with StrictLogging {

  private val drsScheme: String = fileSystemProvider.getScheme

  override def name: String = "DRS"

  override def build(pathAsString: String, pathBuilders: PathBuilders): Try[Path] = {
    if (pathAsString.startsWith(s"$drsScheme://")) {
      Try(createDrsOrOtherPath(pathAsString, pathBuilders))
    } else {
      Failure(new IllegalArgumentException(s"$pathAsString does not have a $drsScheme scheme."))
    }
  }

  private def createDrsOrOtherPath(pathAsString: String, pathBuilders: PathBuilders): Path = {
    def drsPath = DrsPath(fileSystemProvider.getCloudNioPath(pathAsString), requesterPaysProjectIdOption)
    if (preResolve) {
      maybeCreateOtherPath(pathAsString, pathBuilders) getOrElse drsPath
    } else {
      drsPath
    }
  }

  private def maybeCreateOtherPath(pathAsString: String, pathBuilders: PathBuilders): Option[Path] = {

    def logAttempt[A](description: String, attempt: => A): Option[A] = {
      logTry(description, Try(attempt))
    }

    def logTry[A](description: String, tried: Try[A]): Option[A] = {
      tried match {
        case Success(result) => Option(result)
        case Failure(exception) =>
          logFailure(description, exception)
      }
    }

    def logFailure(description: String, reason: Any): None.type = {
      logger.debug(s"Unable to $description, will use a DrsPath to access '$pathAsString': $reason")
      None
    }

    def logSuccess(path: Path): Option[Unit] = {
      logger.debug(s"Converted '$pathAsString' to '$path'")
      Option(())
    }

    val gcsPathOption = for {
      // Right now, only pre-resolving GCS. In the future, could pull others like FTP, HTTP, S3, OSS, SRA, etc.
      gsUriOption <- logAttempt(
        "resolve the uri through Martha",
        DrsResolver.getSimpleGsUri(pathAsString, fileSystemProvider.drsPathResolver).unsafeRunSync()
      )
      gcsUrlWithNoCreds <- gsUriOption
      gcsPath <- logAttempt(
        s"create a GcsPath for '$gcsUrlWithNoCreds'",
        PathFactory.buildPath(gcsUrlWithNoCreds, pathBuilders),
      )
      // Extra: Make sure the GcsPath _actually_ has permission to access the path
      _ <- logAttempt(s"access '$gcsPath' with GCS credentials", gcsPath.size)
      _ <- logSuccess(gcsPath)
    } yield gcsPath

    gcsPathOption
  }
}
