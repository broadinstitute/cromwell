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

import scala.concurrent.ExecutionContext
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

    val (newBatch, resultF) = batch
      .queue(abortQuery)
    resultF
      .map {
        case Success(Right(operation @ _)) =>
          println(s"AbortRequestHandler: operation succeeded ${operation.getName}")
          abortQuery.requester ! BatchAbortRequestSuccessful(abortQuery.jobId.jobId)
          Success(())
        case Success(Left(job @ _)) =>
          println("AbortRequestHandler: operation failed due to no operation object")
          // TODO: Alex - we can likely avoid this by using generics on the callback object
          onFailure(
            new SystemBatchApiException(
              new RuntimeException(
                "This is likely a programming error, onSuccess was called without an Operation object"
              )
            )
          )
        case Failure(ex: BatchApiException) =>
          println(s"AbortRequestHandler: operation failed (BatchApiException) -${ex.getMessage}")
          onFailure(ex)
        case Failure(ex) =>
          println(s"AbortRequestHandler: operation failed (unknown) - ${ex.getMessage}")
          onFailure(new SystemBatchApiException(ex))
      }
      .onComplete {
        case Success(_) =>
        case Failure(ex) =>
          println(s"AbortRequestHandler: operation failed (global) - ${ex.getMessage}")
          onFailure(new SystemBatchApiException(ex))
      }

    newBatch
  }
}
