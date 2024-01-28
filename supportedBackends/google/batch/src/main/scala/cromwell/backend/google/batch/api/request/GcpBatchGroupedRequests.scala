package cromwell.backend.google.batch.api.request

import com.google.cloud.batch.v1.{BatchServiceClient, BatchServiceSettings}
import cromwell.backend.google.batch.api.BatchApiRequestManager._

import scala.util.{Failure, Success}
import scala.util.control.NonFatal

// Like PAPIv2 JsonBatchCallback
// TODO: Alex - this can likely be removed
trait OperationCallback {
  def onSuccess(request: BatchApiRequest, jobName: String, operation: Option[com.google.longrunning.Operation]): Unit
  def onFailure(error: Throwable): Unit
}

class GcpBatchGroupedRequests(batchSettings: BatchServiceSettings) {
  private var _requests: List[(BatchApiRequest, OperationCallback)] = List.empty
  def queue(that: BatchApiRequest, callback: OperationCallback): Unit = {
    println("GcpBatchGroupedRequests: Enqueue request")
    _requests = (that, callback) :: _requests
  }
//  def requests: List[BatchApiRequest] = _requests.reverse.map(_._1)

  def size: Int = _requests.size

  def execute(): Unit = {
    println(s"GcpBatchGroupedRequests: execute ${size} requests")
    // TODO: Alex - handle errors?
    val result = scala.util.Using(BatchServiceClient.create(batchSettings)) { client =>
      // TODO: Alex - Verify whether this already handles retries
      // TODO: Alex - complete the responses with the result
      _requests.map { case (request, callback) =>
        try {
          val (jobName, operation) = request match {
            case r: BatchStatusPollRequest =>
              client.getJob(r.httpRequest).getName -> None

            case r: BatchRunCreationRequest =>
              client.createJob(r.httpRequest).getName -> None

            case r: BatchAbortRequest =>
              r.httpRequest.getName -> Option(client.deleteJobCallable().call(r.httpRequest))
          }
          callback.onSuccess(request, jobName, operation)
        } catch {
          case NonFatal(ex) => callback.onFailure(ex)
        }
      }
    }

    result match {
      // TODO: Alex - revisit this
      case Failure(exception) =>
        println(s"GcpBatchGroupedRequests: failed to execute ${size} requests: ${exception.getMessage}")
        throw new RuntimeException("Failed", exception)
      case Success(_) =>
        println(s"GcpBatchGroupedRequests: successfully executed ${size} requests")
        ()
    }
  }
}
