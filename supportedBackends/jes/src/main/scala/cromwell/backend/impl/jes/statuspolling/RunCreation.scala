package cromwell.backend.impl.jes.statuspolling

import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.{GoogleJsonError, GoogleJsonErrorContainer}
import com.google.api.client.http.{HttpHeaders, HttpRequest}
import com.google.api.services.genomics.model.Operation
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager._
import cromwell.backend.standard.StandardAsyncJob

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

private[statuspolling] trait RunCreation { this: JesPollingActor =>

  private def runCreationResultHandler(originalRequest: JesApiQuery, completionPromise: Promise[Try[Unit]]) = new JsonBatchCallback[Operation] {
    override def onSuccess(operation: Operation, responseHeaders: HttpHeaders): Unit = {
      originalRequest.requester ! getJob(operation)
      completionPromise.trySuccess(Success(()))
      ()
    }

    override def onFailure(e: GoogleJsonError, responseHeaders: HttpHeaders): Unit = {
      pollingManager ! JesApiRunCreationQueryFailed(originalRequest, new JesApiException(GoogleJsonException(e, responseHeaders)))
      completionPromise.trySuccess(Failure(new Exception(mkErrorString(e))))
      ()
    }
  }

  def enqueueRunCreationInBatch(runCreationQuery: JesRunCreationQuery, batch: BatchRequest): Future[Try[Unit]] = {
    val completionPromise = Promise[Try[Unit]]()
    val resultHandler = runCreationResultHandler(runCreationQuery, completionPromise)
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
