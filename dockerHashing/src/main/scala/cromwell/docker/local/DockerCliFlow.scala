package cromwell.docker.local

import java.util.concurrent.TimeoutException

import cats.effect.{ContextShift, IO, Timer}
import cromwell.docker.DockerInfoActor._
import cromwell.docker._
import org.http4s.client.Client

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.util.{Failure, Success}

/**
  * A docker flow using the CLI to return docker hashes.
  */
class DockerCliFlow(implicit ec: ExecutionContext) extends DockerRegistry {
  implicit val cs = IO.contextShift(ec)
  implicit val timer = IO.timer(ec)

  // If the docker cli hangs it would be difficult to debug. So timeout the first request after a short duration.
  // https://github.com/docker/docker/issues/18279
  // https://github.com/docker/docker/issues/12606
  lazy val firstLookupTimeout = 5.seconds

  override def accepts(dockerImageIdentifier: DockerImageIdentifier): Boolean = true

  override def run(dockerInfoContext: DockerInfoContext)(implicit client: Client[IO]) = {
    implicit val timer = IO.timer(ec)

    DockerCliFlow.lookupHashOrTimeout(firstLookupTimeout)(dockerInfoContext)
      .flatMap({
        // If the image isn't there, pull it and try again
        case (_: DockerInfoNotFound, _) =>
          DockerCliFlow.pull(dockerInfoContext)
          DockerCliFlow.lookupHashOrTimeout(firstLookupTimeout)(dockerInfoContext)
        case other => IO.pure(other)
      })
  }

  override def config = DockerRegistryConfig.default
}

object DockerCliFlow {
  /**
    * Lookup the hash for the image referenced in the context.
    *
    * @param context The image to lookup.
    * @return The docker hash response plus the context of our flow.
    */
  private def lookupHash(context: DockerInfoContext): (DockerInfoResponse, DockerInfoContext) = {
    val dockerCliKey = cliKeyFromImageId(context)
    DockerInfoActor.logger.debug("Looking up hash of {}", dockerCliKey.fullName)
    val result = DockerCliClient.lookupHash(dockerCliKey) match {
      case Success(None) => DockerInfoNotFound(context.request)
      case Success(Some(hash)) => DockerHashResult.fromString(hash) match {
        case Success(r) => DockerInfoSuccessResponse(DockerInformation(r, None), context.request)
        case Failure(t) => DockerInfoFailedResponse(t, context.request)
      }
      case Failure(throwable) => DockerInfoFailedResponse(throwable, context.request)
    }
    // give the compiler a hint on the debug() override we're trying to use.
    DockerInfoActor.logger.debug("Hash result of {} was {}", dockerCliKey.fullName, result.asInstanceOf[Any])
    (result, context)
  }

  /**
    * Lookup the hash for the image referenced in the context within the timeout, or DockerHashFailedResponse
    * containing a TimeoutException.
    *
    * @param timeout   How long to wait for the exception.
    * @param context   The image to lookup.
    * @return The docker hash response plus the context of our flow.
    */
  private def lookupHashOrTimeout(timeout: FiniteDuration)
                                 (context: DockerInfoContext)
                                 (implicit cs: ContextShift[IO], timer: Timer[IO]): IO[(DockerInfoResponse, DockerInfoContext)] = {
    IO(lookupHash(context)).timeout(timeout)
      .handleErrorWith({
        case _: TimeoutException => IO.pure {
          val dockerCliKey = cliKeyFromImageId(context)
          val exception = new TimeoutException(
            s"""|Timeout while looking up hash of ${dockerCliKey.fullName}.
                |Ensure that docker is running correctly.
                |""".stripMargin)
          DockerInfoFailedResponse(exception, context.request) -> context
        }
        case other => IO.pure(DockerInfoFailedResponse(other, context.request) -> context)
      })
  }

  /**
    * Pull the docker image referenced in context.
    *
    * @param context The image to pull.
    * @return The context of our flow.
    */
  def pull(context: DockerInfoContext): DockerInfoContext = {
    val dockerCliKey = cliKeyFromImageId(context)
    DockerInfoActor.logger.info(s"Attempting to pull {}", dockerCliKey.fullName)
    val result = DockerCliClient.pull(dockerCliKey)
    result match {
      case Success(_) => DockerInfoActor.logger.info("Pulled {}", dockerCliKey.fullName)
      case Failure(throwable) => DockerInfoActor.logger.error("Docker pull failed", throwable)
    }
    context
  }

  /** Utility for converting the flow image id to the format output by the docker cli. */
  private def cliKeyFromImageId(context: DockerInfoContext): DockerCliKey = {
    val imageId = context.dockerImageID
    (imageId.host, imageId.repository) match {
      case (None, None) =>
        // For docker hub images (host == None), and don't include "library".
        val repository = imageId.image
        val tag = imageId.reference
        DockerCliKey(repository, tag)
      case _ =>
        // For all other images, include the host and repository.
        val repository = s"${imageId.hostAsString}${imageId.nameWithDefaultRepository}"
        val tag = imageId.reference
        DockerCliKey(repository, tag)
    }
  }
}
