package cromwell.backend.google.batch.api.request

import akka.actor.ActorRef
import com.google.cloud.batch.v1.CreateJobRequest
import cromwell.backend.google.batch.api.BatchApiRequestManager._
import cromwell.backend.standard.StandardAsyncJob

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

trait RunRequestHandler { this: RequestHandler =>
  private def runCreationResultHandler(originalRequest: BatchApiRequest,
                                       completionPromise: Promise[Try[Unit]],
                                       pollingManager: ActorRef
  ) = new OperationCallback {
    override def onSuccess(operation: BatchApiRequest): Unit = {
      originalRequest.requester ! getJob(operation)
      completionPromise.trySuccess(Success(()))
      ()
    }

    override def onFailure(e: String): Unit = {
      // TODO: Alex - find a better way to report errors
      val rootCause = new Exception(e)
//      val rootCause = new Exception(mkErrorString(e))

      // TODO: Alex - differentiate between system and user errors
      // See com.google.api.gax.rpc.ApiException
      val failureException = new SystemBatchApiException(rootCause)
//      val failureException = if (e.getCode.toString.startsWith(HttpUserErrorCodeInitialNumber)) {
//        val helpfulHint = if (rootCause.getMessage.contains("unsupported accelerator")) {
//          Option("See https://cloud.google.com/compute/docs/gpus/ for a list of supported accelerators.")
//        } else None
//
//        new UserPAPIApiException(GoogleJsonException(e, responseHeaders), helpfulHint)
//      } else {
//        new SystemPAPIApiException(GoogleJsonException(e, responseHeaders))
//      }

      pollingManager ! BatchApiRunCreationQueryFailed(originalRequest, failureException)
      completionPromise.trySuccess(Failure(failureException))
      ()
    }
  }

  def handleRequest(runCreationQuery: BatchRunCreationRequest,
                    batch: GcpBatchGroupedRequests,
                    pollingManager: ActorRef
  ): Future[Try[Unit]] = {
    val completionPromise = Promise[Try[Unit]]()
    val resultHandler = runCreationResultHandler(runCreationQuery, completionPromise, pollingManager)
    addRunCreationToBatch(runCreationQuery.httpRequest, batch, resultHandler)
    completionPromise.future
  }

  private def addRunCreationToBatch(request: CreateJobRequest,
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

  // TODO: Alex - replace operation type by the batch one
  private def getJob(operation: BatchApiRequest) = StandardAsyncJob(operation.getName)
}

object RunRequestHandler {
  // Because HTTP 4xx errors indicate user error:
  val HttpUserErrorCodeInitialNumber: String = "4"
}
