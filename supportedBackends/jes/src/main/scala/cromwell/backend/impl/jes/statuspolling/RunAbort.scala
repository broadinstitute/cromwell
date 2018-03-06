package cromwell.backend.impl.jes.statuspolling

import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.{GoogleJsonError, GoogleJsonErrorContainer}
import com.google.api.client.http.{HttpHeaders, HttpRequest}
import com.google.api.services.genomics.model.Operation
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager._
import cromwell.backend.impl.jes.statuspolling.RunAbort.JesAbortRequestSuccessful

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

private[statuspolling] trait RunAbort extends PapiInstrumentation { this: JesPollingActor =>

  private def abortResultHandler(originalRequest: JesAbortQuery, completionPromise: Promise[Try[Unit]]) = new JsonBatchCallback[Operation] {
    override def onSuccess(operation: Operation, responseHeaders: HttpHeaders): Unit = {
      abortSuccess()
      originalRequest.requester ! getJob(operation)
      completionPromise.trySuccess(Success(()))
      ()
    }

    override def onFailure(e: GoogleJsonError, responseHeaders: HttpHeaders): Unit = {
      // No need to fail the request if the job was already cancelled, we're all good
      if (Option(e.getCode).contains(400) && Option(e.getMessage).contains("Operation has already been canceled")) {
        originalRequest.requester ! JesAbortRequestSuccessful(originalRequest.run.job.jobId)
        completionPromise.trySuccess(Success(()))
      } else {
        pollingManager ! JesApiAbortQueryFailed(originalRequest, new JesApiException(GoogleJsonException(e, responseHeaders)))
        completionPromise.trySuccess(Failure(new Exception(mkErrorString(e))))
      }

      ()
    }
  }

  def enqueueAbortInBatch(abortQuery: JesAbortQuery, batch: BatchRequest): Future[Try[Unit]] = {
    val completionPromise = Promise[Try[Unit]]()
    val resultHandler = abortResultHandler(abortQuery, completionPromise)
    addAbortToBatch(abortQuery.httpRequest, batch, resultHandler)
    completionPromise.future
  }

  private def addAbortToBatch(request: HttpRequest, batch: BatchRequest, resultHandler: JsonBatchCallback[Operation]): Unit = {
    /*
      * Manually enqueue the request instead of doing it through the RunPipelineRequest
      * as it would unnecessarily rebuild the request (which we already have)
     */
    batch.queue(request, classOf[Operation], classOf[GoogleJsonErrorContainer], resultHandler)
    ()
  }

  private def getJob(operation: Operation) = JesAbortRequestSuccessful(operation.getName)
}

object RunAbort {
  case class JesAbortRequestSuccessful(operationId: String)
}
