package cromwell.backend.impl.jes.statuspolling

import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpHeaders
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model.{Operation, RunPipelineRequest}
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
      pollingManager ! JesApiQueryCreationFailed(originalRequest, new JesApiException(GoogleJsonException(e, responseHeaders)))
      completionPromise.trySuccess(Failure(new Exception(mkErrorString(e))))
      ()
    }
  }

  def enqueueRunCreationInBatch(runCreationQuery: JesRunCreationQuery, batch: BatchRequest): Future[Try[Unit]] = {
    val completionPromise = Promise[Try[Unit]]()
    val resultHandler = runCreationResultHandler(runCreationQuery, completionPromise)
    addRunCreationToBatch(runCreationQuery.rpr, runCreationQuery.genomicsInterface, batch, resultHandler)
    completionPromise.future
  }

  private def addRunCreationToBatch(rpr: RunPipelineRequest, genomicsInterface: Genomics, batch: BatchRequest, resultHandler: JsonBatchCallback[Operation]) =
    genomicsInterface.pipelines().run(rpr).queue(batch, resultHandler)

  private def getJob(operation: Operation) = StandardAsyncJob(operation.getName)
}
