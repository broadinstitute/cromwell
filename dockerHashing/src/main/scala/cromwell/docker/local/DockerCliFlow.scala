package cromwell.docker.local

import java.util.concurrent.TimeoutException

import akka.actor.Scheduler
import akka.stream.scaladsl.{Flow, GraphDSL, Merge, Partition}
import akka.stream.{ActorMaterializer, FlowShape}
import cromwell.docker.DockerHashActor._
import cromwell.docker.{DockerFlow, DockerHashActor, DockerHashResult, DockerImageIdentifierWithoutHash}

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

/**
  * A docker flow using the CLI to return docker hashes.
  */
class DockerCliFlow(implicit ec: ExecutionContext, materializer: ActorMaterializer, scheduler: Scheduler)
  extends DockerFlow {

  // If the docker cli hangs it would be difficult to debug. So timeout the first request after a short duration.
  // https://github.com/docker/docker/issues/18279
  // https://github.com/docker/docker/issues/12606
  lazy val firstLookupTimeout = 5.seconds

  override def accepts(dockerImageIdentifierWithoutHash: DockerImageIdentifierWithoutHash): Boolean = true

  override def buildFlow() = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._

    // Try a lookup.
    val firstLookup = builder.add(
      Flow
        .fromFunction(DockerCliFlow.lookupHashOrTimeout(firstLookupTimeout))
        .mapAsync(1)(identity)
    )

    // Check if the hash was not found
    val partition = builder.add(Partition[(DockerHashResponse, DockerHashContext)](2, {
      case (_: DockerHashNotFound, _) => 0 // If we didn't find the docker hash, go to out(0)
      case _ => 1 // Otherwise, go to out(1)
    }))

    // If the hash wasn't found, `docker pull` then try a lookup a second time. If the docker pull fails for any
    // reason, the error is ignored, allowing the second lookup to return not found.
    val dockerPull = builder.add(Flow.fromFunction(Function.tupled(DockerCliFlow.pull _)))
    val secondLookup = builder.add(Flow.fromFunction(DockerCliFlow.lookupHash))

    // Use either the first or second docker lookup.
    val merge = builder.add(Merge[(DockerHashResponse, DockerHashContext)](2))

    // @formatter:off
    firstLookup.out ~> partition
                       partition.out(0) ~> dockerPull ~> secondLookup ~> merge.in(0)
                       partition.out(1)               ~>                 merge.in(1)
    // @formatter:on

    FlowShape(firstLookup.in, merge.out)
  }
}

object DockerCliFlow {
  /**
    * Lookup the hash for the image referenced in the context.
    *
    * @param context The image to lookup.
    * @return The docker hash response plus the context of our flow.
    */
  private def lookupHash(context: DockerHashContext): (DockerHashResponse, DockerHashContext) = {
    val dockerCliKey = cliKeyFromImageId(context)
    DockerHashActor.logger.debug("Looking up hash of {}", dockerCliKey.fullName)
    val result = DockerCliClient.lookupHash(dockerCliKey) match {
      case Success(None) => DockerHashNotFound(context.request)
      case Success(Some(hash)) => DockerHashResponseSuccess(DockerHashResult(hash), context.request)
      case Failure(throwable) => DockerHashFailedResponse(throwable, context.request)
    }
    // give the compiler a hint on the debug() override we're trying to use.
    DockerHashActor.logger.debug("Hash result of {} was {}", dockerCliKey.fullName, result.asInstanceOf[Any])
    (result, context)
  }

  /**
    * Lookup the hash for the image referenced in the context within the timeout, or DockerHashFailedResponse
    * containing a TimeoutException.
    *
    * @param timeout   How long to wait for the exception.
    * @param context   The image to lookup.
    * @param scheduler Schedules the timout exception.
    * @return The docker hash response plus the context of our flow.
    */
  private def lookupHashOrTimeout(timeout: FiniteDuration)
                                 (context: DockerHashContext)
                                 (implicit ec: ExecutionContext,
                                  scheduler: Scheduler): Future[(DockerHashResponse, DockerHashContext)] = {
    val normal = Future(lookupHash(context))
    val delayed = akka.pattern.after(timeout, scheduler) {
      val dockerCliKey = cliKeyFromImageId(context)
      val exception = new TimeoutException(
        s"""|Timeout while looking up hash of ${dockerCliKey.fullName}.
            |Ensure that docker is running correctly.
            |""".stripMargin)
      val response = DockerHashFailedResponse(exception, context.request)
      Future.successful((response, context))
    }
    Future.firstCompletedOf(Seq(normal, delayed))
  }

  /**
    * Pull the docker image referenced in context.
    *
    * @param notUsed Output of `lookupHash`, passed in only to help creating the graph easier.
    * @param context The image to pull.
    * @return The context of our flow.
    */
  private def pull(notUsed: DockerHashResponse, context: DockerHashContext): DockerHashContext = {
    val dockerCliKey = cliKeyFromImageId(context)
    DockerHashActor.logger.info(s"Attempting to pull {}", dockerCliKey.fullName)
    val result = DockerCliClient.pull(dockerCliKey)
    result match {
      case Success(_) => DockerHashActor.logger.info("Pulled {}", dockerCliKey.fullName)
      case Failure(throwable) => DockerHashActor.logger.error("Docker pull failed", throwable)
    }
    context
  }

  /** Utility for converting the flow image id to the format output by the docker cli. */
  private def cliKeyFromImageId(context: DockerHashContext): DockerCliKey = {
    val imageId = context.dockerImageID
    (imageId.host, imageId.repository) match {
      case (None, "library") =>
        // For docker hub images (host == None), and don't include "library".
        val repository = imageId.image
        val tag = imageId.reference
        DockerCliKey(repository, tag)
      case _ =>
        // For all other images, include the host and repository.
        val repository = s"${imageId.hostAsString}${imageId.repository}/${imageId.image}"
        val tag = imageId.reference
        DockerCliKey(repository, tag)
    }
  }
}
