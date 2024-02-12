package cromwell.backend.google.batch.api.request

import akka.actor.ActorRef
import com.google.cloud.batch.v1.Job
import com.google.longrunning.Operation
import com.typesafe.scalalogging.LazyLogging
import cromwell.backend.google.batch.actors.BatchApiAbortClient.BatchAbortRequestSuccessful
import cromwell.backend.google.batch.api.BatchApiRequestManager.{
  BatchAbortRequest,
  BatchApiAbortQueryFailed,
  BatchApiException,
  BatchApiRequest,
  SystemBatchApiException
}

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

trait AbortRequestHandler extends LazyLogging { this: RequestHandler =>
  private def abortRequestResultHandler(
    originalRequest: BatchAbortRequest,
    completionPromise: Promise[Try[Unit]],
    pollingManager: ActorRef
  ) = new OperationCallback {
    override def onSuccess(request: BatchApiRequest, result: Either[Job, Operation]): Unit = {
      result match {
        case Right(operation @ _) =>
          originalRequest.requester ! BatchAbortRequestSuccessful(originalRequest.jobId.jobId)
          completionPromise.trySuccess(Success(()))

        case Left(job @ _) =>
          // TODO: Alex - we can likely avoid this by using generics on the callback object
          onFailure(
            new SystemBatchApiException(
              new RuntimeException(
                "This is likely a programming error, onSuccess was called without an Operation object"
              )
            )
          )
      }
      ()
    }

    override def onFailure(ex: BatchApiException): Unit = {
      pollingManager ! BatchApiAbortQueryFailed(originalRequest, ex)
      completionPromise.trySuccess(Failure(ex))
      ()
    }
  }

  def handleRequest(
    abortQuery: BatchAbortRequest,
    batch: GcpBatchGroupedRequests,
    pollingManager: ActorRef
  ): Future[Try[Unit]] = {
    val completionPromise = Promise[Try[Unit]]()
    val resultHandler = abortRequestResultHandler(abortQuery, completionPromise, pollingManager)
    addAbortQueryToBatch(abortQuery, batch, resultHandler)
    completionPromise.future
  }

  private def addAbortQueryToBatch(
    request: BatchAbortRequest,
    batch: GcpBatchGroupedRequests,
    resultHandler: OperationCallback
  ): Unit = {
    /*
     * Manually enqueue the request instead of doing it through the RunPipelineRequest
     * as it would unnecessarily rebuild the request (which we already have)
     */
    batch.queue(request, resultHandler)
    ()
  }
}
