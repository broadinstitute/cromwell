package cromwell.filesystems.drs

import java.net.URI

import cats.data.NonEmptyList
import cloud.nio.impl.drs.DrsCloudNioFileSystemProvider
import com.google.common.net.UrlEscapers
import cromwell.core.path.PathFactory.PathBuilders
import cromwell.core.path.{Path, PathBuilder, PathFactory}
import org.slf4j.Logger

import scala.util.{Failure, Success, Try}

case class DrsPathBuilder(fileSystemProvider: DrsCloudNioFileSystemProvider,
                          requesterPaysProjectIdOption: Option[String]) extends PathBuilder {

  private val drsScheme: String = fileSystemProvider.getScheme
  private val GcsScheme: String = "gs"

  override def name: String = "DRS"

  /**
    * Unlike other cloud providers where directories are identified with a trailing slash at the end like `gs://bucket/dir/`,
    * DRS has a concept of bundles for directories (not supported yet). This sanitizes the path by removing the trailing '/'(s)
    * so that the CloudNioFileSystemProvider & CloudNioPath does not treat a DRS path ending with '/' as a directory, otherwise
    * its size is returned as 0
    */
  private def removeTrailingSlashes(path: String): String = {
    if (path.length > (drsScheme.length + 3)) { //3: length of '://'
      val pathArray = path.split(s"://")
      val transformedPath = pathArray(1).replaceAll("[/]+$", "")
      s"$drsScheme://$transformedPath"
    }
    else path
  }

  override def build(logger: Logger, pathAsString: String, pathBuilders: PathBuilders): Try[Path] = {
    if (pathAsString.startsWith(s"$drsScheme://")) {
      val pathWithoutTrailingSlashes = removeTrailingSlashes(pathAsString)
      Try(URI.create(UrlEscapers.urlFragmentEscaper().escape(pathWithoutTrailingSlashes))) flatMap { uri =>
        if (!Option(uri.getScheme).exists(_.equalsIgnoreCase(fileSystemProvider.getScheme))) {
          Failure(new IllegalArgumentException(s"$pathAsString does not have a $drsScheme scheme."))
        } else {
          // Right now, only pre-resolving GCS. In the future, could pull others like HTTP, AWS, etc.
          val gcsPathBuilders = pathBuilders.filter(_.name == "Google Cloud Storage")
          Try(createDrsOrOtherPath(logger, uri, gcsPathBuilders))
        }
      }
    } else {
      Failure(new IllegalArgumentException(s"$pathAsString does not have a $drsScheme scheme."))
    }
  }

  private def createDrsOrOtherPath(logger: Logger, uri: URI, pathBuilders: PathBuilders): Path = {

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
      logger.info(s"Unable to $description, will use a DrsPath to access '$uri': $reason")
      None
    }

    def logSuccess(path: Path): Option[Unit] = {
      logger.info(s"Converted '$uri' to '$path'")
      Option(())
    }

    val gcsPathOption = for {
      // If passed an empty list of path builders, returns a None and does not run the conversion
      nonEmptyPathBuilders <- NonEmptyList.fromList(pathBuilders)
      accessToken <- logAttempt("get an access token", fileSystemProvider.fileProvider.getAccessToken())
      marthaResponse <- logAttempt(
        "resolve the uri through Martha",
        fileSystemProvider.drsPathResolver.resolveDrsThroughMartha(uri.toString, Option(accessToken)).unsafeRunSync(),
      )
      // Right now, only pre-resolving GCS. In the future, could pull others like FTP, HTTP, S3, OSS, SRA, etc.
      gcsUrlWithNoCreds <- marthaResponse.googleServiceAccount match {
        case Some(_) => logFailure("convert uri", "requires additional auth")
        case None => logAttempt(
          s"retrieve $GcsScheme url from martha response",
          marthaResponse.gsUri.get,
        )
      }
      gcsPath <- logAttempt(
        s"create a GcsPath for '$gcsUrlWithNoCreds'",
        PathFactory.buildPath(logger, gcsUrlWithNoCreds, nonEmptyPathBuilders.toList),
      )
      // Make sure the GcsPath _actually_ has permission to access the path
      _ <- logAttempt(s"access '$gcsPath' with GCS credentials", gcsPath.size)
      _ <- logSuccess(gcsPath)
    } yield gcsPath

    gcsPathOption getOrElse DrsPath(fileSystemProvider.getPath(uri), requesterPaysProjectIdOption)
  }
}
