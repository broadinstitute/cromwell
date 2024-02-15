package cromwell.backend.google.batch.api.request

import akka.actor.ActorRef
import com.typesafe.scalalogging.LazyLogging
import cromwell.backend.google.batch.actors.BatchApiAbortClient.BatchAbortRequestSuccessful
import cromwell.backend.google.batch.api.BatchApiRequestManager.{
  BatchAbortRequest,
  BatchApiAbortQueryFailed,
  BatchApiException,
  SystemBatchApiException
}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.control.NonFatal
import scala.util.{Failure, Success, Try}

trait AbortRequestHandler extends LazyLogging { this: RequestHandler =>

  def handleRequest(
    abortQuery: BatchAbortRequest,
    batch: GcpBatchGroupedRequests,
    pollingManager: ActorRef
  )(implicit ec: ExecutionContext): Future[Try[Unit]] = {
    def onFailure(ex: BatchApiException): Try[Unit] = {
      pollingManager ! BatchApiAbortQueryFailed(abortQuery, ex)
      Failure(ex)
    }

    batch
      .queue(abortQuery)
      .map {
        case Success(Right(operation @ _)) =>
          abortQuery.requester ! BatchAbortRequestSuccessful(abortQuery.jobId.jobId)
          Success(())
        case Success(Left(job @ _)) =>
          // TODO: Alex - we can likely avoid this by using generics on the callback object
          onFailure(
            new SystemBatchApiException(
              new RuntimeException(
                "This is likely a programming error, onSuccess was called without an Operation object"
              )
            )
          )
        case Failure(ex: BatchApiException) => onFailure(ex)
        case Failure(ex) => onFailure(new SystemBatchApiException(ex))
      }
      .recover { case NonFatal(ex) =>
        onFailure(new SystemBatchApiException(ex))
      }
  }
}
