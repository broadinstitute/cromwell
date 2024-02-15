package cromwell.backend.google.batch.api.request

import akka.actor.ActorRef
import cromwell.backend.google.batch.api.BatchApiRequestManager.{
  BatchApiException,
  BatchApiStatusQueryFailed,
  BatchStatusPollRequest,
  SystemBatchApiException
}

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success, Try}

trait GetRequestHandler { this: RequestHandler =>

  def handleRequest(pollRequest: BatchStatusPollRequest, batch: GcpBatchGroupedRequests, pollingManager: ActorRef)(
    implicit ec: ExecutionContext
  ): GcpBatchGroupedRequests = {
    def onFailure(ex: BatchApiException): Try[Unit] = {
      pollingManager ! BatchApiStatusQueryFailed(pollRequest, ex)
      Failure(ex)
    }

    println("GetRequestHandler: enqueue request")
    val (newBatch, resultF) = batch
      .queue(pollRequest)

    resultF
      .map {
        case Success(Right(operation @ _)) =>
          println("GetRequestHandler: operation failed due to no job object")
          onFailure(
            new SystemBatchApiException(
              new RuntimeException(
                "This is likely a programming error, onSuccess was called without a job object"
              )
            )
          )

        case Success(Left(job @ _)) =>
          println(s"GetRequestHandler: operation succeeded ${job.getName}")
          pollRequest.requester ! job
          Success(())

        case Failure(ex: BatchApiException) =>
          println(s"GetRequestHandler: operation failed (BatchApiException) -${ex.getMessage}")
          onFailure(ex)
        case Failure(ex) =>
          println(s"GetRequestHandler: operation failed (unknown) - ${ex.getMessage}")
          onFailure(new SystemBatchApiException(ex))
      }
      .onComplete {
        case Success(_) =>
        case Failure(ex) =>
          println(s"GetRequestHandler: operation failed (global) - ${ex.getMessage}")
          onFailure(new SystemBatchApiException(ex))
      }

    newBatch
  }
}
