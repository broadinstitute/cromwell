package cromwell.backend.google.batch.api.request

import akka.actor.ActorRef
import com.typesafe.scalalogging.LazyLogging
import cromwell.backend.google.batch.actors.BatchApiAbortClient
import cromwell.backend.google.batch.actors.BatchApiAbortClient.BatchAbortRequestSuccessful
import cromwell.backend.google.batch.api.BatchApiRequestManager.{
  BatchAbortRequest,
  BatchApiAbortQueryFailed,
  BatchApiException,
  SystemBatchApiException
}
import cromwell.backend.google.batch.api.BatchApiResponse

import scala.concurrent.ExecutionContext
import scala.util.control.NonFatal
import scala.util.{Failure, Success, Try}

trait AbortRequestHandler extends LazyLogging { this: RequestHandler =>

  def handleRequest(
    abortQuery: BatchAbortRequest,
    batch: GcpBatchGroupedRequests,
    pollingManager: ActorRef
  )(implicit ec: ExecutionContext): GcpBatchGroupedRequests = {
    def onFailure(ex: BatchApiException): Try[Unit] = {
      pollingManager ! BatchApiAbortQueryFailed(abortQuery, ex)
      Failure(ex)
    }

    val (newBatch, resultF) = batch.enqueue(abortQuery)

    val _ = resultF
      .map {
        case Success(BatchApiResponse.DeleteJobRequested(BatchApiAbortClient.BatchAbortRequestSuccessful(jobId))) =>
          logger.info(s"Abort job requested successfully: $jobId")
          abortQuery.requester ! BatchAbortRequestSuccessful(jobId)
          Success(())

        case Success(BatchApiResponse.DeleteJobRequested(BatchApiAbortClient.BatchOperationIsAlreadyTerminal(jobId))) =>
          logger.info(s"Job was not aborted because it is already in a terminal state: $jobId")
          abortQuery.requester ! BatchAbortRequestSuccessful(jobId)
          Success(())

        case Success(result) =>
          logger.error(
            s"Programming error, abort operation failed due to no DeleteJobRequested object, got this instead: $result"
          )
          onFailure(
            new SystemBatchApiException(
              new RuntimeException(
                "This is likely a programming error, onSuccess was called without a DeleteJobRequested object"
              )
            )
          )

        case Failure(ex: BatchApiException) =>
          logger.error(s"Abort operation failed: ${abortQuery.jobId.jobId}", ex)
          onFailure(ex)

        case Failure(ex) =>
          logger.error(s"Abort operation failed (unknown reason): ${abortQuery.jobId.jobId}", ex)
          onFailure(new SystemBatchApiException(ex))
      }
      .recover { case NonFatal(ex) =>
        logger.error(s"Abort operation failed (global handler): ${abortQuery.jobId.jobId}", ex)
        onFailure(new SystemBatchApiException(ex))
      }

    newBatch
  }
}
