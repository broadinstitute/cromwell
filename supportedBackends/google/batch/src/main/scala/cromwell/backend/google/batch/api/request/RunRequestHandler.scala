package cromwell.backend.google.batch.api.request

import akka.actor.ActorRef
import com.typesafe.scalalogging.LazyLogging
import cromwell.backend.google.batch.api.BatchApiRequestManager._
import cromwell.backend.google.batch.api.BatchApiResponse
import cromwell.backend.standard.StandardAsyncJob

import scala.concurrent.ExecutionContext
import scala.util.control.NonFatal
import scala.util.{Failure, Success, Try}

trait RunRequestHandler extends LazyLogging { this: RequestHandler =>

  def handleRequest(runCreationQuery: BatchRunCreationRequest,
                    batch: GcpBatchGroupedRequests,
                    pollingManager: ActorRef
  )(implicit ec: ExecutionContext): GcpBatchGroupedRequests = {
    def onFailure(ex: BatchApiException): Try[Unit] = {
      pollingManager ! BatchApiRunCreationQueryFailed(runCreationQuery, ex)
      Failure(ex)
    }

    val (newBatch, resultF) = batch.queue(runCreationQuery)
    val _ = resultF
      .map {
        case Success(BatchApiResponse.JobCreated(job)) =>
          logger.info(s"Run operation succeeded ${job.getName}")
          runCreationQuery.requester ! getJob(job.getName)
          Success(())

        case Success(result) =>
          logger.error(s"Run operation failed due to no job object. got this instead: ${result}")
          onFailure(
            new SystemBatchApiException(
              new RuntimeException(
                "This is likely a programming error, onSuccess was called without a job object"
              )
            )
          )

        case Failure(ex: BatchApiException) =>
          logger.error(s"Run operation failed (BatchApiException)", ex)
          onFailure(ex)

        case Failure(ex) =>
          logger.error(s"Run operation failed (unknown)", ex)
          onFailure(new SystemBatchApiException(ex))
      }
      .recover { case NonFatal(ex) =>
        logger.error(s"Run operation failed (global)", ex)
        onFailure(new SystemBatchApiException(ex))
      }

    newBatch
  }

  private def getJob(jobName: String) = StandardAsyncJob(jobName)
}
