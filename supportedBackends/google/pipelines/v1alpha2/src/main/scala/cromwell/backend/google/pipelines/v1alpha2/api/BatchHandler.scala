package cromwell.backend.google.pipelines.v1alpha2.api

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpRequest
import com.google.api.services.genomics.Genomics
import cromwell.backend.google.pipelines.common.api.PipelinesApiBatchHandler
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.{PAPIAbortRequest, PAPIApiRequest, PAPIRunCreationRequest, PAPIStatusPollRequest}
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode

import scala.collection.JavaConverters._

class BatchHandler(applicationName: String) extends PipelinesApiBatchHandler 
  with CreateRequestBatchHandler
  with StatusRequestBatchHandler 
  with AbortRequestBatchHandler {

  override def makeBatchRequest = new Genomics.Builder(GoogleAuthMode.httpTransport, GoogleAuthMode.jsonFactory, (_: HttpRequest) => ())
    .setApplicationName(applicationName)
    .build().batch()

  override def enqueue[T <: PAPIApiRequest](papiApiRequest: T, batchRequest: BatchRequest, pollingManager: ActorRef) = papiApiRequest match {
    case create: PAPIRunCreationRequest => enqueueRunCreationInBatch(create, batchRequest, pollingManager)
    case status: PAPIStatusPollRequest => enqueueStatusPollInBatch(status, batchRequest, pollingManager)
    case abort: PAPIAbortRequest => enqueueAbortInBatch(abort, batchRequest, pollingManager)
  }
  
  protected def mkErrorString(e: GoogleJsonError) = e.getErrors.asScala.toList.mkString(", ")
}
