package cromwell.backend.google.pipelines.v2alpha1.api.request

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.{GoogleJsonError, GoogleJsonErrorContainer}
import com.google.api.client.http.{HttpHeaders, HttpRequest}
import com.google.api.services.genomics.v2alpha1.model.Operation
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager._
import cromwell.backend.standard.StandardAsyncJob
import cromwell.backend.google.pipelines.v2alpha1.api.request.RunRequestHandler.HttpUserErrorCodeInitialNumber

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

trait RunRequestHandler { this: RequestHandler =>
  private def runCreationResultHandler(originalRequest: PAPIApiRequest, completionPromise: Promise[Try[Unit]], pollingManager: ActorRef) = new JsonBatchCallback[Operation] {
    override def onSuccess(operation: Operation, responseHeaders: HttpHeaders): Unit = {
      originalRequest.requester ! getJob(operation)
      completionPromise.trySuccess(Success(()))
      ()
    }

    override def onFailure(e: GoogleJsonError, responseHeaders: HttpHeaders): Unit = {

      val rootCause = new Exception(mkErrorString(e))

      val papiFailureException = if (e.getCode.toString.startsWith(HttpUserErrorCodeInitialNumber)) {
        val helpfulHint = if (rootCause.getMessage.contains("unsupported accelerator")) {
          Option("See https://cloud.google.com/compute/docs/gpus/ for a list of supported accelerators.")
        } else None

        new UserPAPIApiException(GoogleJsonException(e, responseHeaders), helpfulHint)
      } else {
        new SystemPAPIApiException(GoogleJsonException(e, responseHeaders))
      }

      pollingManager ! PipelinesApiRunCreationQueryFailed(originalRequest, papiFailureException)
      completionPromise.trySuccess(Failure(papiFailureException))
      ()
    }
  }

  def handleRequest(runCreationQuery: PAPIRunCreationRequest, batch: BatchRequest, pollingManager: ActorRef): Future[Try[Unit]] = {
    val completionPromise = Promise[Try[Unit]]()
    val resultHandler = runCreationResultHandler(runCreationQuery, completionPromise, pollingManager)
    addRunCreationToBatch(runCreationQuery.httpRequest, batch, resultHandler)
    completionPromise.future
  }

  private def addRunCreationToBatch(request: HttpRequest, batch: BatchRequest, resultHandler: JsonBatchCallback[Operation]): Unit = {
    /*
      * Manually enqueue the request instead of doing it through the RunPipelineRequest
      * as it would unnecessarily rebuild the request (which we already have)
     */
    batch.queue(request, classOf[Operation], classOf[GoogleJsonErrorContainer], resultHandler)
    ()
  }

  private def getJob(operation: Operation) = StandardAsyncJob(operation.getName)
}

object RunRequestHandler {
  // Because HTTP 4xx errors indicate user error:
  val HttpUserErrorCodeInitialNumber: String = "4"
}
