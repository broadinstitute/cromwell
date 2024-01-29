package cromwell.backend.google.batch.api.request

import akka.actor.ActorRef
//import com.google.api.client.googleapis.batch.BatchRequest
//import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.gax.rpc.FixedHeaderProvider
import com.google.cloud.batch.v1.BatchServiceSettings
import com.google.common.collect.ImmutableMap
import cromwell.backend.google.batch.api.BatchApiRequestManager
import cromwell.backend.google.batch.api.BatchApiRequestManager._
import cromwell.backend.google.batch.models.GcpBatchConfigurationAttributes.BatchRequestTimeoutConfiguration
//import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import org.slf4j.{Logger, LoggerFactory}

//import java.net.URL
//import scala.jdk.CollectionConverters._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

trait BatchApiRequestHandler {
  def makeBatchRequest: GcpBatchGroupedRequests

  def enqueue[T <: BatchApiRequestManager.BatchApiRequest](request: T,
                                                           batchRequest: GcpBatchGroupedRequests,
                                                           pollingManager: ActorRef
  )(implicit ec: ExecutionContext): Future[Try[Unit]]
}

object RequestHandler {
  val logger: Logger = LoggerFactory.getLogger("BatchApiRequestHandler")
}

class RequestHandler(applicationName: String,
//                     endpointUrl: URL,
                     batchRequestTimeoutConfiguration: BatchRequestTimeoutConfiguration
) extends BatchApiRequestHandler
    with RunRequestHandler
    with GetRequestHandler
    with AbortRequestHandler {

  override def makeBatchRequest: GcpBatchGroupedRequests = {
    val headers = ImmutableMap.of("user-agent", "cromwell")
    val headerProvider = FixedHeaderProvider.create(headers)
    val batchSettings = BatchServiceSettings.newBuilder.setHeaderProvider(headerProvider).build
    new GcpBatchGroupedRequests(batchSettings)
  }

  override def enqueue[T <: BatchApiRequestManager.BatchApiRequest](request: T,
                                                                    batchRequest: GcpBatchGroupedRequests,
                                                                    pollingManager: ActorRef
  )(implicit ec: ExecutionContext): Future[Try[Unit]] = request match {
    case create: BatchRunCreationRequest =>
      // TODO: Alex - Remove this
      println(ec.hashCode())
      handleRequest(create, batchRequest, pollingManager)
    case status: BatchStatusPollRequest =>
      handleRequest(status, batchRequest, pollingManager)

    case abort: BatchAbortRequest =>
      handleRequest(abort, batchRequest, pollingManager)
  }

//  protected def mkErrorString(e: GoogleJsonError): String = e.getErrors.asScala.toList.mkString(", ")
}
