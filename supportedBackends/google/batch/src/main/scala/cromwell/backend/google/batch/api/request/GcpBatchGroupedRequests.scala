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
class GcpBatchGroupedRequests(batchSettings: BatchServiceSettings,
                              requests: List[(BatchApiRequest, Promise[Try[Either[Job, Operation]]])] = List.empty
) {

//  private var _requests: List[(BatchApiRequest, Promise[Try[Either[Job, Operation]]])] = List.empty
  def queue(that: BatchApiRequest): (GcpBatchGroupedRequests, Future[Try[Either[Job, Operation]]]) = {
    // TODO: Alex - remove log
    println(s"GcpBatchGroupedRequests: Enqueue request, queue has ${size} items")

    val promise = Promise[Try[Either[Job, Operation]]]()
    val newRequests = (that, promise) :: requests
    new GcpBatchGroupedRequests(batchSettings, newRequests) -> promise.future
  }
//  def requests: List[BatchApiRequest] = _requests.reverse.map(_._1)

  def size: Int = requests.size

  private def internalExecute(client: BatchServiceClient, request: BatchApiRequest): Try[Either[Job, Operation]] =
    try {
      val result = request match {
        case r: BatchStatusPollRequest =>
          // TODO: Alex - remove log
          println(s"GcpBatchGroupedRequests: getJob")
          Left(client.getJob(r.httpRequest))

        case r: BatchRunCreationRequest =>
          // TODO: Alex - remove log
          println(s"GcpBatchGroupedRequests: createJob")
          Left(client.createJob(r.httpRequest))

        case r: BatchAbortRequest =>
          // TODO: Alex - we should find a way to detect when the operation is already terminal
          // then, throw BatchOperationIsAlreadyTerminal
          // different to PAPIv2, this call does not abort the job but deletes it, so, we need to be careful
          // to not delete jobs in a terminal state

          // TODO: Alex - remove log
          println(s"GcpBatchGroupedRequests: deleteJob")

          Right(client.deleteJobCallable().call(r.httpRequest))
      }
      println(s"GcpBatchGroupedRequests: result = ${result.map(_.getName).left.map(_.getName)}")
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
        println(s"GcpBatchGroupedRequests: failure = ${apiException.getMessage}")
        Failure(failureException)

      case NonFatal(ex) =>
        println(s"GcpBatchGroupedRequests: failure (nonfatal) = ${ex.getMessage}")
        Failure(new SystemBatchApiException(ex))
    }

  // TODO: Alex - add retries and the logic from BatchRequest
  def execute(implicit ec: ExecutionContext): Future[List[Try[Unit]]] = {
    println(s"GcpBatchGroupedRequests: execute ${size} requests")

    val client = BatchServiceClient.create(batchSettings)

    // TODO: Alex - Verify whether this already handles retries
    // TODO: The sdk calls seem to be blocking, consider using a Future + scala.concurrent.blocking
    //       Check whether sending many requests in parallel could be an issue
    val futures = requests.reverse.map { case (request, promise) =>
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
      println("GcpBatchGroupedRequests: closing client")
      try
        client.close()
      catch {
        case NonFatal(ex) =>
          println(s"GcpBatchGroupedRequests: failed to close client${ex.getMessage}")
      }
    }

    result
  }
}
