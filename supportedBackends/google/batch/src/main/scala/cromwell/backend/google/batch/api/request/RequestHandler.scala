package cromwell.backend.google.batch.api.request

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.gax.rpc.FixedHeaderProvider
import com.google.api.services.lifesciences.v2beta.CloudLifeSciences
import com.google.cloud.batch.v1.BatchServiceSettings
import com.google.common.collect.ImmutableMap
import cromwell.backend.google.batch.api.BatchApiRequestManager
import cromwell.backend.google.batch.api.BatchApiRequestManager._
import cromwell.backend.google.batch.models.GcpBatchConfigurationAttributes.BatchRequestTimeoutConfiguration
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import org.slf4j.{Logger, LoggerFactory}

import java.net.URL
import scala.jdk.CollectionConverters._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

object RequestHandler {
  val logger: Logger = LoggerFactory.getLogger("BatchApiRequestHandler")
}

class RequestHandler(applicationName: String,
                     endpointUrl: URL,
                     batchRequestTimeoutConfiguration: BatchRequestTimeoutConfiguration
) /*extends PipelinesApiRequestHandler with*/
    extends RunRequestHandler
    with GetRequestHandler
    with AbortRequestHandler {

  def makeBatchRequest: GcpBatchGroupedRequests = {
    val headers = ImmutableMap.of("user-agent", "cromwell")
    val headerProvider = FixedHeaderProvider.create(headers)
    val batchSettings = BatchServiceSettings.newBuilder.setHeaderProvider(headerProvider).build
    new GcpBatchGroupedRequests(batchSettings)
  }

  def enqueue[T <: BatchApiRequestManager.BatchApiRequest](papiApiRequest: T,
                                                           batchRequest: GcpBatchGroupedRequests,
                                                           pollingManager: ActorRef
  )(implicit ec: ExecutionContext): Future[Try[Unit]] = papiApiRequest match {
    case create: BatchRunCreationRequest => handleRequest(create, batchRequest, pollingManager)
    case status: BatchStatusPollRequest => handleRequest(status, batchRequest, pollingManager)
    case abort: BatchAbortRequest => handleRequest(abort, batchRequest, pollingManager)
  }

  private def withClient[T](f: BatchServiceClient => T): T = {
    // set user agent to cromwell so requests can be differentiated on batch
    val headers = ImmutableMap.of("user-agent", "cromwell")
    val headerProvider = FixedHeaderProvider.create(headers)
    val batchSettings = BatchServiceSettings.newBuilder.setHeaderProvider(headerProvider).build
    val client = BatchServiceClient.create(batchSettings)
    try
      f(client)
    finally
      client.close()
  }
//  protected def mkErrorString(e: GoogleJsonError): String = e.getErrors.asScala.toList.mkString(", ")
}
