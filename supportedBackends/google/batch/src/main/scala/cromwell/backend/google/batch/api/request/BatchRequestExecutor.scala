package cromwell.backend.google.batch.api.request

import com.google.api.gax.rpc.{ApiException, StatusCode}
import com.google.cloud.batch.v1._
import com.typesafe.scalalogging.LazyLogging
import cromwell.backend.google.batch.actors.BatchApiAbortClient.{
  BatchAbortRequestSuccessful,
  BatchOperationIsAlreadyTerminal
}
import cromwell.backend.google.batch.api.BatchApiRequestManager._
import cromwell.backend.google.batch.api.{BatchApiRequestManager, BatchApiResponse}
import cromwell.backend.google.batch.models.RunStatus
import cromwell.core.ExecutionEvent
import io.grpc.Status
import org.apache.commons.lang3.exception.ExceptionUtils

import scala.annotation.unused
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.jdk.CollectionConverters.ListHasAsScala
import scala.util.control.NonFatal
import scala.util.{Failure, Success, Try}

trait BatchRequestExecutor {
  def execute(groupedRequests: GcpBatchGroupedRequests)(implicit ec: ExecutionContext): Future[List[Try[Unit]]]
}

object BatchRequestExecutor {

  class CloudImpl(batchSettings: BatchServiceSettings) extends BatchRequestExecutor with LazyLogging {

    // TODO: Alex - add retries and the logic from BatchRequest
    def execute(groupedRequests: GcpBatchGroupedRequests)(implicit ec: ExecutionContext): Future[List[Try[Unit]]] = {
      val requests = groupedRequests.entries
      logger.info(s"Execute ${requests.size} requests")

      if (requests.isEmpty) Future.successful(List.empty)
      else nonEmptyExecute(requests)
    }

    private def nonEmptyExecute(
      requests: List[(BatchApiRequestManager.BatchApiRequest, Promise[Try[BatchApiResponse]])]
    )(implicit ec: ExecutionContext): Future[List[Try[Unit]]] = {
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
        try
          client.close()
        catch {
          case NonFatal(ex) =>
            logger.warn(s"Failed to close batch client: ${ex.getMessage}", ex)
        }
      }

      result
    }

    private def internalExecute(client: BatchServiceClient, request: BatchApiRequest): Try[BatchApiResponse] =
      try
        request match {
          case r: BatchStatusPollRequest =>
            val result = internalGetHandler(client, r.httpRequest)
            Success(result)

          case r: BatchRunCreationRequest =>
            val result = BatchApiResponse.JobCreated(client.createJob(r.httpRequest))
            Success(result)

          case r: BatchAbortRequest =>
            // different to PAPIv2, this call does not abort the job but deletes it, so, we need to be careful
            // to not delete jobs in a terminal state
            val getResult = internalGetHandler(client, GetJobRequest.newBuilder.setName(r.httpRequest.getName).build())
            val abortResult = getResult.status match {
              case _: RunStatus.TerminalRunStatus => BatchOperationIsAlreadyTerminal(r.jobId.jobId)
              case _ =>
                // After playing with the sdk, it seems that operation.getResultCase is always RESULT_NOT_SET
                // When there is an error, an exception is thrown right away.
                @unused
                val operation = client.deleteJobCallable().call(r.httpRequest)
                BatchAbortRequestSuccessful(r.jobId.jobId)
            }

            Success(BatchApiResponse.DeleteJobRequested(abortResult))
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
          Failure(failureException)

        case NonFatal(ex) => Failure(new SystemBatchApiException(ex))
      }

    // A different handler is used when fetching the job status to map exceptions to the correct RunStatus
    private def internalGetHandler(client: BatchServiceClient, request: GetJobRequest): BatchApiResponse.StatusQueried =
      try {
        val job = client.getJob(request)
        val result = interpretOperationStatus(job)
        BatchApiResponse.StatusQueried(result)
      } catch {
        // We don't need to detect preemptible VMs because that's handled automatically by GCP
        case apiException: ApiException if apiException.getStatusCode.getCode == StatusCode.Code.RESOURCE_EXHAUSTED =>
          BatchApiResponse.StatusQueried(RunStatus.AwaitingCloudQuota)
      }

    private[request] def interpretOperationStatus(job: Job): RunStatus =
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
      /* TODO: This is an example printing the events from GCP, do we need to do anything else in the mapping?
Event type=STATUS_CHANGED
time=seconds: 1712173852,nanos: 952604950
taskState=STATE_UNSPECIFIED,
description=Job state is set from QUEUED to SCHEDULED for job projects/392615380452/locations/us-south1/jobs/job-ba81bad8-82e9-4d95-8fc0-04dfbbd746da.
taskExecution.exitCode=0

Event type=STATUS_CHANGED,
time=seconds: 1712173947, nanos: 568998105
taskState=STATE_UNSPECIFIED
description=Job state is set from SCHEDULED to RUNNING for job projects/392615380452/locations/us-south1/jobs/job-ba81bad8-82e9-4d95-8fc0-04dfbbd746da.
taskExecution.exitCode=0

Event type=STATUS_CHANGED
time=seconds: 1712173989, nanos: 937816549
taskState=STATE_UNSPECIFIED
description=Job state is set from RUNNING to SUCCEEDED for job projects/392615380452/locations/us-south1/jobs/job-ba81bad8-82e9-4d95-8fc0-04dfbbd746da.
taskExecution.exitCode=0
       */
      events.map { e =>
        /* TODO: This is an example output, verify that this is what we need
ExecutionEvent(Job state is set from QUEUED to SCHEDULED for job projects/392615380452/locations/us-south1/jobs/job-321db1bc-9a68-4171-aa2a-46885d781656.,2024-04-03T20:10:01.704137839Z,None)
Event - ExecutionEvent(Job state is set from SCHEDULED to RUNNING for job projects/392615380452/locations/us-south1/jobs/job-321db1bc-9a68-4171-aa2a-46885d781656.,2024-04-03T20:11:30.631264449Z,None)
Event - ExecutionEvent(Job state is set from RUNNING to SUCCEEDED for job projects/392615380452/locations/us-south1/jobs/job-321db1bc-9a68-4171-aa2a-46885d781656.,2024-04-03T20:12:16.898798407Z,None)
         */
        val time = java.time.Instant
          .ofEpochSecond(e.getEventTime.getSeconds, e.getEventTime.getNanos.toLong)
          .atOffset(java.time.ZoneOffset.UTC)
        ExecutionEvent(name = e.getDescription, offsetDateTime = time)
      }
  }
}
