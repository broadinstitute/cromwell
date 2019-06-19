package cromwell.backend.google.pipelines.v1alpha2.api.request

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.{GoogleJsonError, GoogleJsonErrorContainer}
import com.google.api.client.http.{HttpHeaders, HttpRequest}
import com.google.api.services.genomics.model.Operation
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager._
import cromwell.backend.google.pipelines.common.api.clients.PipelinesApiAbortClient.{PAPIAbortRequestSuccessful, PAPIOperationAlreadyCancelled, PAPIOperationHasAlreadyFinished}

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

trait AbortRequestHandler { this: RequestHandler =>
  private def abortResultHandler(originalRequest: PAPIAbortRequest, completionPromise: Promise[Try[Unit]], pollingManager: ActorRef) = new JsonBatchCallback[Operation] {
    override def onSuccess(operation: Operation, responseHeaders: HttpHeaders): Unit = {
      originalRequest.requester ! PAPIAbortRequestSuccessful(originalRequest.jobId.jobId)
      completionPromise.trySuccess(Success(()))
      ()
    }

    override def onFailure(e: GoogleJsonError, responseHeaders: HttpHeaders): Unit = {
      // No need to fail the request if the job was already cancelled, we're all good
      if (Option(e.getCode).contains(400) && Option(e.getMessage).contains("Operation has already been canceled")) {
        originalRequest.requester ! PAPIOperationAlreadyCancelled(originalRequest.jobId.jobId)
        completionPromise.trySuccess(Success(()))
      } else if (Option(e.getCode).contains(400) && Option(e.getMessage).contains("Operation has already finished")) {
        originalRequest.requester ! PAPIOperationHasAlreadyFinished(originalRequest.jobId.jobId)
        completionPromise.trySuccess(Success(()))
      } else {
        pollingManager ! PipelinesApiAbortQueryFailed(originalRequest, new SystemPAPIApiException(GoogleJsonException(e, responseHeaders)))
        completionPromise.trySuccess(Failure(new Exception(mkErrorString(e))))
      }

      ()
    }
  }

  def enqueueAbortInBatch(abortQuery: PAPIAbortRequest, batch: BatchRequest, pollingMangaer: ActorRef): Future[Try[Unit]] = {
    val completionPromise = Promise[Try[Unit]]()
    val resultHandler = abortResultHandler(abortQuery, completionPromise, pollingMangaer)
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
}
