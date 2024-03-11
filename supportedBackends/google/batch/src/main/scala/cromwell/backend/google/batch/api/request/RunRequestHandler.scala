package cromwell.backend.google.batch.api.request

import akka.actor.ActorRef
import cromwell.backend.google.batch.api.BatchApiRequestManager._
import cromwell.backend.standard.StandardAsyncJob

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success, Try}

trait RunRequestHandler { this: RequestHandler =>

  def handleRequest(runCreationQuery: BatchRunCreationRequest,
                    batch: GcpBatchGroupedRequests,
                    pollingManager: ActorRef
  )(implicit ec: ExecutionContext): GcpBatchGroupedRequests = {
    def onFailure(ex: BatchApiException): Try[Unit] = {
      pollingManager ! BatchApiRunCreationQueryFailed(runCreationQuery, ex)
      Failure(ex)
    }

    println("RunRequestHandler: enqueue request")
    val (newBatch, resultF) = batch.queue(runCreationQuery)
    resultF
      .map {
        case Success(BatchResponse.JobCreated(job)) =>
          println(s"RunRequestHandler: operation succeeded ${job.getName}")
          runCreationQuery.requester ! getJob(job.getName)
          Success(())

        case Success(_) =>
          println("RunRequestHandler: operation failed due to no job object")
          onFailure(
            new SystemBatchApiException(
              new RuntimeException(
                "This is likely a programming error, onSuccess was called without a job object"
              )
            )
          )

        case Failure(ex: BatchApiException) =>
          println(s"RunRequestHandler: operation failed (BatchApiException) -${ex.getMessage}")
          onFailure(ex)

        case Failure(ex) =>
          println(s"RunRequestHandler: operation failed (unknown) - ${ex.getMessage}")
          onFailure(new SystemBatchApiException(ex))
      }
      .onComplete {
        case Success(_) =>
        case Failure(ex) =>
          println(s"RunRequestHandler: operation failed (global) - ${ex.getMessage}")
          onFailure(new SystemBatchApiException(ex))
      }
    newBatch
  }

  private def getJob(jobName: String) = StandardAsyncJob(jobName)
}
