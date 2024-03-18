package cromwell.backend.google.batch.api.request

import akka.actor.ActorRef
import com.google.api.gax.rpc.FixedHeaderProvider
import com.google.cloud.batch.v1.BatchServiceSettings
import com.google.common.collect.ImmutableMap
import cromwell.backend.google.batch.api.BatchApiRequestManager
import cromwell.backend.google.batch.api.BatchApiRequestManager._

import scala.concurrent.ExecutionContext

trait BatchApiRequestHandler {
  def makeBatchRequest: GcpBatchGroupedRequests

  def enqueue[T <: BatchApiRequestManager.BatchApiRequest](request: T,
                                                           batchRequest: GcpBatchGroupedRequests,
                                                           pollingManager: ActorRef
  )(implicit ec: ExecutionContext): GcpBatchGroupedRequests
}

class RequestHandler
    extends BatchApiRequestHandler
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
  )(implicit ec: ExecutionContext): GcpBatchGroupedRequests = request match {
    case create: BatchRunCreationRequest =>
      handleRequest(create, batchRequest, pollingManager)

    case status: BatchStatusPollRequest =>
      handleRequest(status, batchRequest, pollingManager)

    case abort: BatchAbortRequest =>
      handleRequest(abort, batchRequest, pollingManager)
  }
}
