package cromwell.backend.google.pipelines.v1alpha2.api.request

import java.net.URL

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpRequest
import com.google.api.services.genomics.Genomics
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestHandler
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.{PAPIAbortRequest, PAPIApiRequest, PAPIRunCreationRequest, PAPIStatusPollRequest}
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode

import scala.collection.JavaConverters._
import scala.concurrent.ExecutionContext

class RequestHandler(applicationName: String, endpointUrl: URL) extends PipelinesApiRequestHandler 
  with RunRequestHandler
  with GetRequestHandler 
  with AbortRequestHandler {

  override def makeBatchRequest = new Genomics.Builder(GoogleAuthMode.httpTransport, GoogleAuthMode.jsonFactory, (_: HttpRequest) => ())
    .setApplicationName(applicationName)
    .setRootUrl(endpointUrl.toString)
    .build().batch()

  override def enqueue[T <: PAPIApiRequest](papiApiRequest: T, batchRequest: BatchRequest, pollingManager: ActorRef)
                                           (implicit ec: ExecutionContext)= papiApiRequest match {
    case create: PAPIRunCreationRequest => enqueueRunCreationInBatch(create, batchRequest, pollingManager)
    case status: PAPIStatusPollRequest => enqueueStatusPollInBatch(status, batchRequest, pollingManager)
    case abort: PAPIAbortRequest => enqueueAbortInBatch(abort, batchRequest, pollingManager)
  }
  
  protected def mkErrorString(e: GoogleJsonError) = e.getErrors.asScala.toList.mkString(", ")
}
