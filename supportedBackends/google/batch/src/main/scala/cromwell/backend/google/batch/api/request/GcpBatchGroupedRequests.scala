package cromwell.backend.google.batch.api.request

import com.google.api.gax.rpc.ApiException
import com.google.cloud.batch.v1.{BatchServiceClient, BatchServiceSettings, Job}
import cromwell.backend.google.batch.api.BatchApiRequestManager._

import scala.util.{Failure, Success}
import scala.util.control.NonFatal
import com.google.longrunning.Operation

// Like PAPIv2 JsonBatchCallback
// TODO: Alex - this can likely be removed
trait OperationCallback {
  def onSuccess(request: BatchApiRequest, result: Either[Job, Operation]): Unit
  def onFailure(error: BatchApiException): Unit
}

// Mirrors com.google.api.client.googleapis.batch.BatchRequest
// TODO: Alex - This is not Thread-safe yet
class GcpBatchGroupedRequests(batchSettings: BatchServiceSettings) {

  private var _requests: List[(BatchApiRequest, OperationCallback)] = List.empty
  def queue(that: BatchApiRequest, callback: OperationCallback): Unit = {
    println("GcpBatchGroupedRequests: Enqueue request")
    _requests = (that, callback) :: _requests
  }
//  def requests: List[BatchApiRequest] = _requests.reverse.map(_._1)

  def size: Int = _requests.size

  // TODO: Alex - add retries and the logic from BatchRequest
  def execute(): Unit = {
    println(s"GcpBatchGroupedRequests: execute ${size} requests")
    // TODO: Alex - handle errors?
    val result = scala.util.Using(BatchServiceClient.create(batchSettings)) { client =>
      // TODO: Alex - Verify whether this already handles retries
      // TODO: The sdk calls seem to be blocking, consider using a Future + scala.concurrent.blocking
      //       Check whether sending many requests in parallel could be an issue
      _requests.map { case (request, callback) =>
        try {
          val res = request match {
            case r: BatchStatusPollRequest =>
              Left(client.getJob(r.httpRequest))

            case r: BatchRunCreationRequest =>
              Left(client.createJob(r.httpRequest))

            case r: BatchAbortRequest =>
              Right(client.deleteJobCallable().call(r.httpRequest))
          }
          callback.onSuccess(request, res)
        } catch {
          case apiException: ApiException =>
            // Because HTTP 4xx errors indicate user error:
            val HttpUserErrorCodeInitialNumber: String = "4"
            val failureException =
              if (
                apiException.getStatusCode.getCode.getHttpStatusCode.toString.startsWith(HttpUserErrorCodeInitialNumber)
              ) {
                new UserBatchApiException(GoogleBatchException(apiException), None)
              } else {
                new SystemBatchApiException(GoogleBatchException(apiException))
              }
            callback.onFailure(failureException)

          case NonFatal(ex) => callback.onFailure(new SystemBatchApiException(ex))
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
    // its important to clear the queue after executing the requests
    _requests = List.empty
  }
}
