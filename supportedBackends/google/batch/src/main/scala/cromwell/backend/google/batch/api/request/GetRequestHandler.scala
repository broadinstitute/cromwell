package cromwell.backend.google.batch.api.request

import akka.actor.ActorRef
import com.google.cloud.batch.v1.Job
import com.google.longrunning.Operation
import cromwell.backend.google.batch.api.BatchApiRequestManager.{
  BatchApiException,
  BatchApiRequest,
  BatchApiStatusQueryFailed,
  BatchStatusPollRequest,
  SystemBatchApiException
}

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

trait GetRequestHandler { this: RequestHandler =>
  private def pollRequestResultHandler(originalRequest: BatchApiRequest,
                                       completionPromise: Promise[Try[Unit]],
                                       pollingManager: ActorRef
  ) = new OperationCallback {
    override def onSuccess(request: BatchApiRequest, result: Either[Job, Operation]): Unit = {
      result match {
        case Right(_) =>
          // TODO: Alex - we can likely avoid this by using generics on the callback object
          onFailure(
            new SystemBatchApiException(
              new RuntimeException(
                "This is likely a programming error, onSuccess was called without an the job object"
              )
            )
          )

        case Left(job) =>
          originalRequest.requester ! job
          completionPromise.trySuccess(Success(()))
      }
      ()
    }

    override def onFailure(ex: BatchApiException): Unit = {
      pollingManager ! BatchApiStatusQueryFailed(originalRequest, ex)
      completionPromise.trySuccess(Failure(ex))
      ()
    }
  }

  def handleRequest(pollRequest: BatchStatusPollRequest,
                    batch: GcpBatchGroupedRequests,
                    pollingManager: ActorRef
  ): Future[Try[Unit]] = {
    val completionPromise = Promise[Try[Unit]]()
    val resultHandler = pollRequestResultHandler(pollRequest, completionPromise, pollingManager)
    addPollRequestToBatch(pollRequest, batch, resultHandler)
    completionPromise.future
  }

  private def addPollRequestToBatch(request: BatchStatusPollRequest,
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
