package cromwell.backend.google.batch.api.request

import akka.actor.ActorRef
import cromwell.backend.google.batch.api.BatchApiRequestManager
import cromwell.backend.google.batch.api.BatchApiRequestManager._

import scala.concurrent.ExecutionContext

trait BatchApiRequestHandler {
  def makeBatchRequest: GcpBatchGroupedRequests = GcpBatchGroupedRequests.empty

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
