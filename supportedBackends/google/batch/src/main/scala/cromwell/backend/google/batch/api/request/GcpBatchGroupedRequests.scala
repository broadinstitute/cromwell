package cromwell.backend.google.batch.api.request

import com.google.cloud.batch.v1.{BatchServiceClient, BatchServiceSettings}
import cromwell.backend.google.batch.api.BatchApiRequestManager._

import scala.util.control.NonFatal

// Like PAPIv2 JsonBatchCallback
// TODO: Alex - this can likely be removed
trait OperationCallback {
  def onSuccess(operation: BatchApiRequest): Unit
  def onFailure(error: Throwable): Unit
}

class GcpBatchGroupedRequests(batchSettings: BatchServiceSettings) {
  private var _requests: List[(BatchApiRequest, OperationCallback)] = List.empty
  def queue(that: BatchApiRequest, callback: OperationCallback): Unit = _requests = (that, callback) :: _requests
//  def requests: List[BatchApiRequest] = _requests.reverse.map(_._1)

  def size: Int = _requests.size

  def execute(grouped: GcpBatchGroupedRequests): Unit = {
    // TODO: Alex - handle errors?
    val x = scala.util.Using(BatchServiceClient.create(batchSettings)) { client =>
      // TODO: Alex - Verify whether this already handles retries
      // TODO: Alex - complete the responses with the result
      try {}
      _requests.map { case (request, callback) =>
        try {
          val response = request match {
            case r: BatchStatusPollRequest =>
              client.getJob(r.httpRequest)

            case r: BatchRunCreationRequest =>
              client.createJob(r.httpRequest)

            case r: BatchAbortRequest =>
              val x = client.deleteJobCallable().call(r.httpRequest)
          }
          callback.onSuccess(request, response)
        } catch {
          case NonFatal(ex) => callback.onFailure(ex)
        }
      }
    }
  }
}
