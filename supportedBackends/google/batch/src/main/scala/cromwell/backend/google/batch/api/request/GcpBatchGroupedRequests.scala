package cromwell.backend.google.batch.api.request

import com.google.api.gax.rpc.ApiException
import com.google.cloud.batch.v1.{BatchServiceClient, BatchServiceSettings, Job, JobStatus, StatusEvent}
import cromwell.backend.google.batch.api.BatchApiRequestManager._

import scala.util.{Failure, Success, Try}
import scala.util.control.NonFatal
import com.google.longrunning.Operation
import cromwell.backend.google.batch.models.RunStatus
import cromwell.core.ExecutionEvent
import io.grpc.Status
import org.apache.commons.lang3.exception.ExceptionUtils

import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.jdk.CollectionConverters.ListHasAsScala

sealed trait BatchResponse extends Product with Serializable
object BatchResponse {
  case class JobCreated(job: Job) extends BatchResponse
  case class DeleteJobRequested(operation: Operation) extends BatchResponse
  case class StatusQueried(status: RunStatus) extends BatchResponse
}

// Mirrors com.google.api.client.googleapis.batch.BatchRequest
// TODO: Alex - This is not Thread-safe yet
class GcpBatchGroupedRequests(
  batchSettings: BatchServiceSettings,
  requests: List[(BatchApiRequest, Promise[Try[BatchResponse]])] = List.empty
) {

  def queue(that: BatchApiRequest): (GcpBatchGroupedRequests, Future[Try[BatchResponse]]) = {
    // TODO: Alex - remove log
    println(s"GcpBatchGroupedRequests: Enqueue request, queue has ${size} items")

    val promise = Promise[Try[BatchResponse]]()
    val newRequests = (that, promise) :: requests
    new GcpBatchGroupedRequests(batchSettings, newRequests) -> promise.future
  }
//  def requests: List[BatchApiRequest] = _requests.reverse.map(_._1)

  def size: Int = requests.size

  private def internalGetHandler(client: BatchServiceClient, request: BatchStatusPollRequest): Try[BatchResponse] =
    try {
      val job = client.getJob(request.httpRequest)
      val result = interpretOperationStatus(job, request)
      Success(BatchResponse.StatusQueried(result))
    } catch {
      // TODO: We could detect rate limits/quote exceeded, preemptible vms, etc
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

  private[request] def interpretOperationStatus(job: Job, pollingRequest: BatchStatusPollRequest): RunStatus =
    if (Option(job).isEmpty) {
      // TODO: This applies to PAPIv2 but its unlikely to apply to batch
      //
      // It is possible to receive a null via an HTTP 200 with no response. If that happens, handle it and don't crash.
      // https://github.com/googleapis/google-http-java-client/blob/v1.28.0/google-http-client/src/main/java/com/google/api/client/http/HttpResponse.java#L456-L458
      val errorMessage = "Operation returned as empty"
      RunStatus.UnsuccessfulRunStatus(
        Status.UNKNOWN,
        Option(errorMessage),
        Nil,
        None,
        None,
        None,
        wasPreemptible = false
      )
    } else {
      try {
        lazy val events = job.getStatus.getStatusEventsList.asScala.toList
        if (job.getStatus.getState == JobStatus.State.SUCCEEDED) {
          // TODO: We are likely better by removing the params we don't have
          RunStatus.Success(getEventList(events), None, None, None)
        } else if (isQuotaDelayed(events)) {
          RunStatus.AwaitingCloudQuota
        } else if (job.getStatus.getState == JobStatus.State.FAILED) {
          // TODO: How do we get these values? how do we detect is this was a preemptible vm?
          // Status.OK is hardcoded because the request succeeded, we don't have access to the internal response code
          RunStatus.Failed(Status.OK, None, List.empty, getEventList(events), None, None, None)
        } else if (job.getStatus.getState == JobStatus.State.RUNNING) {
          RunStatus.Running
        } else {
          RunStatus.Initializing
        }
      } catch {
        // TODO: Do we care about this?
        case nullPointerException: NullPointerException =>
          throw new RuntimeException(
            s"Caught NPE while interpreting job result ${job.getName}: " +
              s"${ExceptionUtils.getStackTrace(nullPointerException)}. " +
              s"JSON was $job",
            nullPointerException
          )
      }
    }

  private def getEventList(events: List[StatusEvent]): List[ExecutionEvent] =
    // TODO: How do we map these events to the cromwell event type?
    List.empty

  // TODO: How do we detect the quota delayed event? is it even available?
  private def isQuotaDelayed(events: List[StatusEvent]): Boolean =
    false

  private def internalExecute(client: BatchServiceClient, request: BatchApiRequest): Try[BatchResponse] =
    try
      request match {
        case r: BatchStatusPollRequest =>
          // TODO: Alex - remove log
          println(s"GcpBatchGroupedRequests: getJob")
          internalGetHandler(client, r)

        case r: BatchRunCreationRequest =>
          // TODO: Alex - remove log
          println(s"GcpBatchGroupedRequests: createJob")
          val result = BatchResponse.JobCreated(client.createJob(r.httpRequest))
          Success(result)

        case r: BatchAbortRequest =>
          // TODO: Alex - we should find a way to detect when the operation is already terminal
          // then, throw BatchOperationIsAlreadyTerminal
          // different to PAPIv2, this call does not abort the job but deletes it, so, we need to be careful
          // to not delete jobs in a terminal state

          // TODO: Alex - remove log
          println(s"GcpBatchGroupedRequests: deleteJob")
          val result = BatchResponse.DeleteJobRequested(client.deleteJobCallable().call(r.httpRequest))
          Success(result)
      }
    catch {
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
