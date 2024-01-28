package cromwell.backend.google.batch.api.request

import akka.actor.ActorRef
import com.google.cloud.batch.v1.Job
import com.google.longrunning.Operation
//import com.google.api.client.googleapis.batch.BatchRequest
//import com.google.api.client.googleapis.json.GoogleJsonError
//import com.google.api.client.http.HttpHeaders
import com.typesafe.scalalogging.LazyLogging
import cromwell.backend.google.batch.actors.GcpBatchBackendSingletonActor
import cromwell.backend.google.batch.api.BatchApiRequestManager.{
  BatchAbortRequest,
  BatchApiAbortQueryFailed,
  BatchApiRequest,
  SystemBatchApiException
}

//import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
//import org.apache.commons.lang3.StringUtils

import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.util.{Failure, Success, Try}

trait AbortRequestHandler extends LazyLogging { this: RequestHandler =>
  private def abortRequestResultHandler(originalRequest: BatchAbortRequest,
                                        completionPromise: Promise[Try[Unit]],
                                        pollingManager: ActorRef
  ) = new OperationCallback {
    override def onSuccess(request: BatchApiRequest, result: Either[Job, Operation]): Unit = {
      result match {
        case Right(operation) =>
          originalRequest.requester ! GcpBatchBackendSingletonActor.Event.JobAbortRequestSent(operation)
          completionPromise.trySuccess(Success(()))
        case Left(_) =>
          // TODO: Alex - we can likely avoid this by using generics on the callback object
          onFailure(
            new RuntimeException("This is likely a programming error, onSuccess was called without an Operation object")
          )
      }
      ()
    }

    override def onFailure(ex: Throwable): Unit = {
      // TODO: Alex - find a better way to report errors
      val rootCause = ex
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

      pollingManager ! BatchApiAbortQueryFailed(originalRequest, failureException)
      completionPromise.trySuccess(Failure(failureException))
      ()
    }
  }

  def handleRequest(abortQuery: BatchAbortRequest, batch: GcpBatchGroupedRequests, pollingManager: ActorRef)(implicit
    ec: ExecutionContext
  ): Future[Try[Unit]] = {
    println(ec.hashCode())
    val completionPromise = Promise[Try[Unit]]()
    val resultHandler = abortRequestResultHandler(abortQuery, completionPromise, pollingManager)
    addAbortQueryToBatch(abortQuery, batch, resultHandler)
    completionPromise.future
  }

  private def addAbortQueryToBatch(request: BatchAbortRequest,
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

//  def abortJob(jobName: JobName, backendSingletonActor: ActorRef): Unit =
//    backendSingletonActor ! GcpBatchBackendSingletonActor.Action.AbortJob(jobName)
//
//  def abortActorClientReceive: Actor.Receive = {
//    case GcpBatchBackendSingletonActor.Event.JobAbortRequestSent(job) =>
//      log.info(s"Job aborted on GCP: ${job.getName}")
//      abortSuccess()
//
//    case GcpBatchBackendSingletonActor.Event.ActionFailed(jobName, cause) =>
//      val msg = s"Failed to abort job ($jobName) from GCP"
//      log.error(cause, msg)
//      abortFailed()
//  }
}
