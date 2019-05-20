package cromwell.backend.google.pipelines.v2alpha1.api.request

import java.net.URL

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpRequest
import com.google.api.services.genomics.v2alpha1.Genomics
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.{PAPIAbortRequest, PAPIRunCreationRequest, PAPIStatusPollRequest}
import cromwell.backend.google.pipelines.common.api.{PipelinesApiRequestHandler, PipelinesApiRequestManager}
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.concurrent.ExecutionContext

object RequestHandler {
  val logger = LoggerFactory.getLogger("PipelinesApiRequestHandler")
}

class RequestHandler(applicationName: String, endpointUrl: URL, pipelinesApiAttributes: PipelinesApiConfigurationAttributes) extends PipelinesApiRequestHandler
  with RunRequestHandler
  with GetRequestHandler
  with AbortRequestHandler {

  def initializeHttpRequest(r: HttpRequest): Unit = {
    pipelinesApiAttributes.batchRequestTimeoutConfiguration.readTimeoutMillis.foreach {
      t => r.setReadTimeout(t.value)
    }
    pipelinesApiAttributes.batchRequestTimeoutConfiguration.connectTimeoutMillis.foreach {
      t => r.setConnectTimeout(t.value)
    }
  }

  override def makeBatchRequest = new Genomics.Builder(GoogleAuthMode.httpTransport, GoogleAuthMode.jsonFactory, initializeHttpRequest)
    .setApplicationName(applicationName)
    .setRootUrl(endpointUrl.toString)
    .build()
    .batch(initializeHttpRequest)

  override def enqueue[T <: PipelinesApiRequestManager.PAPIApiRequest](papiApiRequest: T, batchRequest: BatchRequest, pollingManager: ActorRef)
                                                                      (implicit ec: ExecutionContext)= papiApiRequest match {
    case create: PAPIRunCreationRequest => handleRequest(create, batchRequest, pollingManager)
    case status: PAPIStatusPollRequest => handleRequest(status, batchRequest, pollingManager)
    case abort: PAPIAbortRequest => handleRequest(abort, batchRequest, pollingManager)
  }

  protected def mkErrorString(e: GoogleJsonError) = e.getErrors.asScala.toList.mkString(", ")
}
