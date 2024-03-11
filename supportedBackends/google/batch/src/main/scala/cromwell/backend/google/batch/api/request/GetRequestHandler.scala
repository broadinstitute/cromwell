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
        case Success(BatchResponse.StatusQueried(status)) =>
          println(s"GetRequestHandler: operation succeeded")
          pollRequest.requester ! status
          Success(())

        case Success(_) =>
          println("GetRequestHandler: operation failed due to no status object")
          onFailure(
            new SystemBatchApiException(
              new RuntimeException(
                "This is likely a programming error, onSuccess was called without a status object"
              )
            )
          )

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
