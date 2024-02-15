package cromwell.backend.google.batch.api.request

import akka.actor.ActorRef
import cromwell.backend.google.batch.api.BatchApiRequestManager._
import cromwell.backend.standard.StandardAsyncJob

import scala.concurrent.{ExecutionContext, Future}
import scala.util.control.NonFatal
import scala.util.{Failure, Success, Try}

trait RunRequestHandler { this: RequestHandler =>

  def handleRequest(runCreationQuery: BatchRunCreationRequest,
                    batch: GcpBatchGroupedRequests,
                    pollingManager: ActorRef
  )(implicit ec: ExecutionContext): Future[Try[Unit]] = {
    def onFailure(ex: BatchApiException): Try[Unit] = {
      pollingManager ! BatchApiRunCreationQueryFailed(runCreationQuery, ex)
      Failure(ex)
    }

    batch
      .queue(runCreationQuery)
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
          runCreationQuery.requester ! getJob(job.getName)
          Success(())

        case Failure(ex: BatchApiException) => onFailure(ex)
        case Failure(ex) => onFailure(new SystemBatchApiException(ex))
      }
      .recover { case NonFatal(ex) =>
        onFailure(new SystemBatchApiException(ex))
      }
  }

  private def getJob(jobName: String) = StandardAsyncJob(jobName)
}
