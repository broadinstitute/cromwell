package cromwell.backend.google.batch.api.request

import akka.actor.ActorRef
import com.typesafe.scalalogging.LazyLogging
import cromwell.backend.google.batch.api.BatchApiRequestManager.{
  BatchApiException,
  BatchApiStatusQueryFailed,
  BatchStatusPollRequest,
  SystemBatchApiException
}
import cromwell.backend.google.batch.api.BatchApiResponse

import scala.concurrent.ExecutionContext
import scala.util.control.NonFatal
import scala.util.{Failure, Success, Try}

trait GetRequestHandler extends LazyLogging { this: RequestHandler =>

  def handleRequest(pollRequest: BatchStatusPollRequest, batch: GcpBatchGroupedRequests, pollingManager: ActorRef)(
    implicit ec: ExecutionContext
  ): GcpBatchGroupedRequests = {
    def onFailure(ex: BatchApiException): Try[Unit] = {
      pollingManager ! BatchApiStatusQueryFailed(pollRequest, ex)
      Failure(ex)
    }

    val (newBatch, resultF) = batch.enqueue(pollRequest)

    val _ = resultF
      .map {
        case Success(BatchApiResponse.StatusQueried(status)) =>
          pollRequest.requester ! status
          Success(())

        case Success(result) =>
          logger.error(
            s"Programming error, get operation failed due to no status object on job ${pollRequest.jobId.jobId}, got this instead: $result"
          )
          onFailure(
            new SystemBatchApiException(
              new RuntimeException(
                "This is likely a programming error, onSuccess was called without a status object"
              )
            )
          )

        case Failure(ex: BatchApiException) =>
          logger.error(s"Get operation failed for job ${pollRequest.jobId.jobId}", ex)
          onFailure(ex)

        case Failure(ex) =>
          logger.error(s"Get operation failed for job ${pollRequest.jobId.jobId} with an unknown reason", ex)
          onFailure(new SystemBatchApiException(ex))
      }
      .recover { case NonFatal(ex) =>
        logger.error(s"Get operation failed (global) for job ${pollRequest.jobId.jobId}", ex)
        onFailure(new SystemBatchApiException(ex))
      }

    newBatch
  }
}
