package cromwell.backend.google.batch.api.request

import com.google.api.gax.rpc.ApiException
import com.google.cloud.batch.v1.{BatchServiceClient, BatchServiceSettings, Job}
import cromwell.backend.google.batch.api.BatchApiRequestManager._

import scala.util.{Failure, Success, Try}
import scala.util.control.NonFatal
import com.google.longrunning.Operation

import scala.concurrent.{ExecutionContext, Future, Promise}

// Mirrors com.google.api.client.googleapis.batch.BatchRequest
// TODO: Alex - This is not Thread-safe yet
class GcpBatchGroupedRequests(batchSettings: BatchServiceSettings) {

  private var _requests: List[(BatchApiRequest, Promise[Try[Either[Job, Operation]]])] = List.empty
  def queue(that: BatchApiRequest): Future[Try[Either[Job, Operation]]] = {
    // TODO: Alex - remove log
    println("GcpBatchGroupedRequests: Enqueue request")

    val promise = Promise[Try[Either[Job, Operation]]]()
    _requests = (that, promise) :: _requests
    promise.future
  }
//  def requests: List[BatchApiRequest] = _requests.reverse.map(_._1)

  def size: Int = _requests.size

  private def internalExecute(client: BatchServiceClient, request: BatchApiRequest): Try[Either[Job, Operation]] =
    try {
      val result = request match {
        case r: BatchStatusPollRequest =>
          Left(client.getJob(r.httpRequest))

        case r: BatchRunCreationRequest =>
          Left(client.createJob(r.httpRequest))

        case r: BatchAbortRequest =>
          // TODO: Alex - we should find a way to detect when the operation is already terminal
          // then, throw BatchOperationIsAlreadyTerminal
          // different to PAPIv2, this call does not abort the job but deletes it, so, we need to be careful
          // to not delete jobs in a terminal state
          Right(client.deleteJobCallable().call(r.httpRequest))
      }
      Success(result)
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
        Failure(failureException)

      case NonFatal(ex) => Failure(new SystemBatchApiException(ex))
    }

  // TODO: Alex - add retries and the logic from BatchRequest
  def execute(implicit ec: ExecutionContext): Future[List[Try[Unit]]] = {
    println(s"GcpBatchGroupedRequests: execute ${size} requests")

    val client = BatchServiceClient.create(batchSettings)

    // its important to clear the queue to not executing the same requests again
    _requests = List.empty

    // TODO: Alex - Verify whether this already handles retries
    // TODO: The sdk calls seem to be blocking, consider using a Future + scala.concurrent.blocking
    //       Check whether sending many requests in parallel could be an issue
    val futures = _requests.map { case (request, promise) =>
      Future {
        scala.concurrent.blocking {
          val result = internalExecute(client, request)
          promise.success(result)
          result
        }
      }
    }

    val result = Future.sequence(futures.map(_.map(_.map(_ => ()))))

    result.onComplete { _ =>
      client.close()
    }

    result
  }
}
