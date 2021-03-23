package cromwell.backend.google.pipelines.v1alpha2.api.request

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.services.genomics.Genomics
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.BatchRequestTimeoutConfiguration
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestHandler
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.{PAPIAbortRequest, PAPIApiRequest, PAPIRunCreationRequest, PAPIStatusPollRequest}
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode

import java.net.URL
import scala.collection.JavaConverters._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

class RequestHandler(applicationName: String,
                     endpointUrl: URL,
                     batchRequestTimeoutConfiguration: BatchRequestTimeoutConfiguration)
  extends PipelinesApiRequestHandler
  with RunRequestHandler
  with GetRequestHandler 
  with AbortRequestHandler {

  override def makeBatchRequest: BatchRequest = {
    val builder =
      new Genomics.Builder(
        GoogleAuthMode.httpTransport,
        GoogleAuthMode.jsonFactory,
        initializeHttpRequest(batchRequestTimeoutConfiguration) _,
      )
        .setApplicationName(applicationName)
        .setRootUrl(endpointUrl.toString)
    val client = builder.build()
    client.batch(client.getRequestFactory.getInitializer)
  }

  override def enqueue[T <: PAPIApiRequest](papiApiRequest: T, batchRequest: BatchRequest, pollingManager: ActorRef)
                                           (implicit ec: ExecutionContext): Future[Try[Unit]] = papiApiRequest match {
    case create: PAPIRunCreationRequest => enqueueRunCreationInBatch(create, batchRequest, pollingManager)
    case status: PAPIStatusPollRequest => enqueueStatusPollInBatch(status, batchRequest, pollingManager)
    case abort: PAPIAbortRequest => enqueueAbortInBatch(abort, batchRequest, pollingManager)
  }
  
  protected def mkErrorString(e: GoogleJsonError): String = e.getErrors.asScala.toList.mkString(", ")
}
