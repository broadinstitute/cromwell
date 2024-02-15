package cromwell.backend.google.batch.api.request

import akka.actor.ActorRef
import cromwell.backend.google.batch.api.BatchApiRequestManager.{
  BatchApiException,
  BatchApiStatusQueryFailed,
  BatchStatusPollRequest,
  SystemBatchApiException
}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.control.NonFatal
import scala.util.{Failure, Success, Try}

trait GetRequestHandler { this: RequestHandler =>

  def handleRequest(pollRequest: BatchStatusPollRequest, batch: GcpBatchGroupedRequests, pollingManager: ActorRef)(
    implicit ec: ExecutionContext
  ): Future[Try[Unit]] = {
    def onFailure(ex: BatchApiException): Try[Unit] = {
      pollingManager ! BatchApiStatusQueryFailed(pollRequest, ex)
      Failure(ex)
    }

    batch
      .queue(pollRequest)
      .map {
        case Success(Right(operation @ _)) =>
          onFailure(
            new SystemBatchApiException(
              new RuntimeException(
                "This is likely a programming error, onSuccess was called without a job object"
              )
            )
          )

        case Success(Left(job @ _)) =>
          pollRequest.requester ! job
          Success(())

        case Failure(ex: BatchApiException) => onFailure(ex)
        case Failure(ex) => onFailure(new SystemBatchApiException(ex))
      }
      .recover { case NonFatal(ex) =>
        onFailure(new SystemBatchApiException(ex))
      }
  }
}
