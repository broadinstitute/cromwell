package cromwell.backend.google.pipelines.v2alpha1.api.request

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.services.genomics.v2alpha1.Genomics
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.BatchRequestTimeoutConfiguration
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.{PAPIAbortRequest, PAPIRunCreationRequest, PAPIStatusPollRequest}
import cromwell.backend.google.pipelines.common.api.{PipelinesApiRequestHandler, PipelinesApiRequestManager}
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import org.slf4j.{Logger, LoggerFactory}

import java.net.URL
import scala.collection.JavaConverters._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

object RequestHandler {
  val logger: Logger = LoggerFactory.getLogger("PipelinesApiRequestHandler")
}

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

  override def enqueue[T <: PipelinesApiRequestManager.PAPIApiRequest](papiApiRequest: T,
                                                                       batchRequest: BatchRequest,
                                                                       pollingManager: ActorRef)
                                                                      (implicit ec: ExecutionContext)
  : Future[Try[Unit]] = papiApiRequest match {
    case create: PAPIRunCreationRequest => handleRequest(create, batchRequest, pollingManager)
    case status: PAPIStatusPollRequest => handleRequest(status, batchRequest, pollingManager)
    case abort: PAPIAbortRequest => handleRequest(abort, batchRequest, pollingManager)
  }

  protected def mkErrorString(e: GoogleJsonError): String = e.getErrors.asScala.toList.mkString(", ")
}
