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
import cromwell.backend.google.batch.models.{GcpBatchExitCode, RunStatus}
import cromwell.core.ExecutionEvent

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

    def execute(groupedRequests: GcpBatchGroupedRequests)(implicit ec: ExecutionContext): Future[List[Try[Unit]]] = {
      val requests = groupedRequests.entries
      logger.info(s"Execute ${requests.size} requests")

      if (requests.isEmpty) Future.successful(List.empty)
      else nonEmptyExecute(requests)
    }

    private def nonEmptyExecute(
      requests: List[(BatchApiRequestManager.BatchApiRequest, Promise[Try[BatchApiResponse]])]
    )(implicit ec: ExecutionContext): Future[List[Try[Unit]]] =
      for {
        client <- Future.fromTry(Try(BatchServiceClient.create(batchSettings)))

        futures = requests.map { case (request, promise) =>
          Future {
            scala.concurrent.blocking {
              val result = internalExecute(client, request)
              promise.success(result)
              result
            }
          }
        }
        resultF = Future.sequence(futures.map(_.map(_.map(_ => ()))))
        _ = resultF.onComplete { _ =>
          try
            client.close()
          catch {
            case NonFatal(ex) =>
              logger.warn(s"Failed to close batch client: ${ex.getMessage}", ex)
          }
        }
        result <- resultF
      } yield result

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
                //
                // TODO: There is a chance we can monitor this operation with
                // client.getHttpJsonOperationsClient.getOperation(operation.getName)
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

        logger.debug(s"Batch Job Status ${job.getName}: ${job.getStatus}")

        val result = interpretOperationStatus(job)
        BatchApiResponse.StatusQueried(result)
      } catch {
        // A job can't be cancelled but deleted, which is why we consider 404 status as the job being cancelled successfully
        case apiException: ApiException if apiException.getStatusCode.getCode == StatusCode.Code.NOT_FOUND =>
          BatchApiResponse.StatusQueried(RunStatus.Aborted)

        // We don't need to detect preemptible VMs because that's handled automatically by GCP
        case apiException: ApiException if apiException.getStatusCode.getCode == StatusCode.Code.RESOURCE_EXHAUSTED =>
          BatchApiResponse.StatusQueried(RunStatus.AwaitingCloudQuota)
      }

    private[request] def interpretOperationStatus(job: Job): RunStatus = {
      lazy val events = getEventList(
        Option(job)
          .flatMap(e => Option(e.getStatus))
          .flatMap(e => Option(e.getStatusEventsList))
          .map(_.asScala.toList)
          .getOrElse(List.empty)
      )
      lazy val exitCode = findBatchExitCode(events)

      if (job.getStatus.getState == JobStatus.State.SUCCEEDED) {
        RunStatus.Success(events)
      } else if (job.getStatus.getState == JobStatus.State.RUNNING) {
        RunStatus.Running
      } else if (job.getStatus.getState == JobStatus.State.FAILED) {
        RunStatus.Failed(exitCode, events)
      } else {
        RunStatus.Initializing
      }
    }

    private def findBatchExitCode(events: List[ExecutionEvent]): Option[GcpBatchExitCode] =
      events.flatMap { e =>
        GcpBatchExitCode.fromEventMessage(e.name.toLowerCase)
      }.headOption

    private def getEventList(events: List[StatusEvent]): List[ExecutionEvent] =
      events.map { e =>
        val time = java.time.Instant
          .ofEpochSecond(e.getEventTime.getSeconds, e.getEventTime.getNanos.toLong)
          .atOffset(java.time.ZoneOffset.UTC)
        ExecutionEvent(name = e.getDescription, offsetDateTime = time)
      }
  }
}
