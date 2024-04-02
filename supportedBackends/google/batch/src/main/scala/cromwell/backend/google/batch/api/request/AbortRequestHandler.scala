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

    val (newBatch, resultF) = batch.queue(abortQuery)

    val _ = resultF
      .map {
        case Success(BatchApiResponse.DeleteJobRequested(operation)) =>
          // After playing with the sdk, it seems that operation.getResultCase is always RESULT_NOT_SET
          // When there is an error, an exception is thrown right away.
          logger.info(s"Operation succeeded ${operation.getName}, result = ${operation.getResultCase}")
          abortQuery.requester ! BatchAbortRequestSuccessful(abortQuery.jobId.jobId)
          Success(())

        case Success(result) =>
          logger.error(s"Abort operation failed due to no operation object, got this instead: $result")
          onFailure(
            new SystemBatchApiException(
              new RuntimeException(
                "This is likely a programming error, onSuccess was called without an Operation object"
              )
            )
          )

        case Failure(ex: BatchApiException) =>
          // TODO: Do we need to do anything about this error? perhaps we could detect whether the job was already in a final state or already deleted
          logger.error(s"Abort operation failed (BatchApiException)", ex)
          onFailure(ex)

        case Failure(ex) =>
          logger.error(s"Abort operation failed (unknown reason)", ex)
          onFailure(new SystemBatchApiException(ex))
      }
      .recover { case NonFatal(ex) =>
        logger.error(s"Abort operation failed (global)", ex)
        onFailure(new SystemBatchApiException(ex))
      }

    newBatch
  }
}
