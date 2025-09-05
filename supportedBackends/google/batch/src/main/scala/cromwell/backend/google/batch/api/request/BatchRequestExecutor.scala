package cromwell.backend.google.batch.api.request

import com.google.api.gax.rpc.StatusCode.Code
import com.google.api.gax.rpc.{ApiException, StatusCode}
import com.google.cloud.batch.v1.AllocationPolicy.ProvisioningModel
import com.google.cloud.batch.v1._
import com.typesafe.scalalogging.LazyLogging
import cromwell.backend.google.batch.actors.BatchApiAbortClient.{
  BatchAbortRequestSuccessful,
  BatchOperationIsAlreadyBeingAborted,
  BatchOperationIsAlreadyTerminal
}
import cromwell.backend.google.batch.api.BatchApiRequestManager._
import cromwell.backend.google.batch.api.{BatchApiRequestManager, BatchApiResponse}
import cromwell.backend.google.batch.models.{GcpBatchExitCode, RunStatus}
import cromwell.core.ExecutionEvent
import cromwell.services.cost.InstantiatedVmInfo
import cromwell.services.metadata.CallMetadataKeys

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
            // In Batch, even if the job is in Cancelled state (one of terminal states), it will still accept
            // a request to cancel the job. So we check the job status before sending the cancel request
            val getResult = internalGetHandler(client, GetJobRequest.newBuilder.setName(r.httpRequest.getName).build())
            val abortResult = getResult.status match {
              // If the job is already in Cancelling state and we try to cancel it again, it will throw a GRPC StatusRuntimeException.
              // Hence in this case we don't send a cancel request again
              case _: RunStatus.Aborting => BatchOperationIsAlreadyBeingAborted(r.jobId.jobId)
              case _: RunStatus.TerminalRunStatus => BatchOperationIsAlreadyTerminal(r.jobId.jobId)
              case _ =>
                // After playing with the sdk, it seems that operation.getResultCase is always RESULT_NOT_SET
                // When there is an error, an exception is thrown right away.
                //
                // TODO: There is a chance we can monitor this operation with
                // client.getHttpJsonOperationsClient.getOperation(operation.getName)
                @unused
                val operation = client.cancelJobCallable().call(r.httpRequest)
                BatchAbortRequestSuccessful(r.jobId.jobId)
            }

            Success(BatchApiResponse.CancelJobRequested(abortResult))
        }
      catch {
        case apiException: ApiException =>
          val exceptionStatusCode = apiException.getStatusCode.getCode
          (exceptionStatusCode, request) match {
            case (Code.ALREADY_EXISTS, req: BatchRunCreationRequest) =>
              val jobName = s"${req.httpRequest.getParent}/jobs/${req.httpRequest.getJobId}"
              logger.info(
                s"Job creation request for $jobName for workflow ${req.workflowId.toString} failed as job already exists. Reconnecting to job instead."
              )
              val result = BatchApiResponse.JobAlreadyExists(jobName)
              Success(result)
            case _ =>
              // Because HTTP 4xx errors indicate user error:
              val HttpUserErrorCodeInitialNumber: String = "4"
              val failureException =
                if (exceptionStatusCode.getHttpStatusCode.toString.startsWith(HttpUserErrorCodeInitialNumber)) {
                  new UserBatchApiException(GoogleBatchException(apiException), None)
                } else {
                  new SystemBatchApiException(GoogleBatchException(apiException))
                }
              Failure(failureException)
          }

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
        // We don't need to detect preemptible VMs because that's handled automatically by GCP
        case apiException: ApiException if apiException.getStatusCode.getCode == StatusCode.Code.RESOURCE_EXHAUSTED =>
          BatchApiResponse.StatusQueried(RunStatus.AwaitingCloudQuota(Seq.empty))
      }

    private[request] def interpretOperationStatus(job: Job): RunStatus = {
      lazy val events = getEventList(
        Option(job)
          .flatMap(e => Option(e.getStatus))
          .flatMap(e => Option(e.getStatusEventsList))
          .map(_.asScala.toList)
          .getOrElse(List.empty)
      )

      // Get vm info for this job
      val allocationPolicy = job.getAllocationPolicy

      // Get instances that can be created with this AllocationPolicy, only instances[0] is supported
      val instancePolicy = allocationPolicy.getInstances(0).getPolicy
      val machineType = instancePolicy.getMachineType

      // SPOT VM is used as preemptible VM instances in Batch. Check for both SPOT or PREEMPTIBLE just to be safe
      val preemptible = (instancePolicy.getProvisioningModelValue == ProvisioningModel.SPOT.getNumber) ||
        (instancePolicy.getProvisioningModelValue == ProvisioningModel.PREEMPTIBLE.getNumber)

      // location list = [regions/us-central1, zones/us-central1-b], region is the first element
      val location = allocationPolicy.getLocation.getAllowedLocationsList.get(0)
      val region =
        if (location.isEmpty)
          "us-central1"
        else
          location.split("/").last

      val instantiatedVmInfo = Some(InstantiatedVmInfo(region, machineType, preemptible))

      job.getStatus.getState match {
        case JobStatus.State.RUNNING => RunStatus.Running(events, instantiatedVmInfo)
        case JobStatus.State.CANCELLATION_IN_PROGRESS => RunStatus.Aborting(events, instantiatedVmInfo)
        case JobStatus.State.SUCCEEDED => RunStatus.Success(events, instantiatedVmInfo)
        case JobStatus.State.FAILED =>
          val batchExitCode =
            events
              .flatMap(e => GcpBatchExitCode.fromEventMessage(e.name))
              .headOption
              .getOrElse(GcpBatchExitCode.Success)
          RunStatus.Failed(batchExitCode, events, instantiatedVmInfo)
        case JobStatus.State.CANCELLED => RunStatus.Aborted(events, instantiatedVmInfo)
        case _ => RunStatus.Initializing(events, instantiatedVmInfo)
      }
    }

    private def getEventList(events: List[StatusEvent]): List[ExecutionEvent] = {
      // on Batch, when job transitions to SCHEDULED state it indicates that the VM is being initialized. Users are billed for this
      // startup time. Hence, the 'vmStartTime' corresponds to when the job enters the SCHEDULED state.
      val startedRegex = ".*to SCHEDULED.*".r

      // job terminal events can occur as below:
      //    - job transitions from a RUNNING state to either SUCCEEDED/FAILED/CANCELLED state
      //    - job never enters RUNNING state and instead transitions from
      //        - SCHEDULED -> FAILED
      //        - SCHEDULED -> SCHEDULED_PENDING_FAILED -> FAILED
      //        - SCHEDULED -> CANCELLATION_IN_PROGRESS -> CANCELLED
      val endedRegex =
        ".*RUNNING to.*|.*SCHEDULED to FAILED.*|.*SCHEDULED_PENDING_FAILED to FAILED.*|.*CANCELLATION_IN_PROGRESS to CANCELLED.*".r

      events.flatMap { e =>
        val time = java.time.Instant
          .ofEpochSecond(e.getEventTime.getSeconds, e.getEventTime.getNanos.toLong)
          .atOffset(java.time.ZoneOffset.UTC)

        // If this event represents the VM starting or stopping, add a special event describing that.
        val startOrEndEvent = (e.getDescription match {
          case startedRegex() => Option(CallMetadataKeys.VmStartTime)
          case endedRegex() => Option(CallMetadataKeys.VmEndTime)
          case _ => None
        }).map(t => ExecutionEvent(name = t, offsetDateTime = time))

        List(ExecutionEvent(name = e.getDescription, offsetDateTime = time)) ++ startOrEndEvent.toList
      }
    }
  }
}
