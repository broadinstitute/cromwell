package cromwell.backend.google.pipelines.v1alpha2.api.request

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.{GoogleJsonError, GoogleJsonErrorContainer}
import com.google.api.client.http.{HttpHeaders, HttpRequest}
import com.google.api.services.genomics.model.Operation
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager._
import cromwell.backend.standard.StandardAsyncJob

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
      pollingManager ! PipelinesApiRunCreationQueryFailed(originalRequest, new SystemPAPIApiException(GoogleJsonException(e, responseHeaders)))
      completionPromise.trySuccess(Failure(new Exception(mkErrorString(e))))
      ()
    }
  }

  def enqueueRunCreationInBatch(runCreationQuery: PAPIRunCreationRequest, batch: BatchRequest, pollingManager: ActorRef): Future[Try[Unit]] = {
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
