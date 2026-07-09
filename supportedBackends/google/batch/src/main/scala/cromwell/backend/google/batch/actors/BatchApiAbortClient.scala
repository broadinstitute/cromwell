package cromwell.backend.google.batch.actors

import akka.actor.{Actor, ActorLogging, ActorRef}
import com.google.cloud.batch.v1.JobName
import cromwell.backend.google.batch.api.BatchApiRequestManager.{BatchAbortRequest, BatchApiAbortQueryFailed}
import cromwell.backend.google.batch.api.GcpBatchRequestFactory
import cromwell.backend.google.batch.monitoring.BatchInstrumentation
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.WorkflowId
import cromwell.core.logging.JobLogging

object BatchApiAbortClient {
  sealed trait BatchAbortRequestSuccess
  case class BatchAbortRequestSuccessful(jobId: String) extends BatchAbortRequestSuccess

  // The operation is no longer running. Maybe it was already cancelled, maybe it finished on its own. We don't know
  // the details and for abort they don't really matter.
  case class BatchOperationIsAlreadyTerminal(jobId: String) extends BatchAbortRequestSuccess

  // The operation is currently being aborted. This means that a previous request to cancel the job was successful,
  // but the job is still in the process of being cancelled.
  case class BatchOperationIsAlreadyBeingAborted(jobId: String) extends BatchAbortRequestSuccess
}

trait BatchApiAbortClient { this: Actor with ActorLogging with JobLogging with BatchInstrumentation =>
  import BatchApiAbortClient._

  def abortJob(workflowId: WorkflowId,
               jobName: JobName,
               backendSingletonActor: ActorRef,
               requestFactory: GcpBatchRequestFactory
  ): Unit =
    backendSingletonActor ! BatchAbortRequest(
      workflowId = workflowId,
      requester = self,
      httpRequest = requestFactory.abortRequest(jobName),
      jobId = StandardAsyncJob(jobName.toString)
    )

  def abortActorClientReceive: Actor.Receive = {
    case BatchAbortRequestSuccessful(jobId) =>
      abortSuccess()
      jobLogger.info(s"Successfully requested cancellation of $jobId")

    // In this case where the job is already being aborted, we log and don't do anything
    case BatchOperationIsAlreadyBeingAborted(jobId) =>
      jobLogger.info(s"Job $jobId is already being cancelled")

    // In this case we could immediately return an aborted handle and spare ourselves a round of polling
    case BatchOperationIsAlreadyTerminal(jobId) =>
      jobLogger.info(s"Job $jobId has already finished")

    case BatchApiAbortQueryFailed(query, e) =>
      jobLogger.error(s"Could not request cancellation of job ${query.jobId}", e)
  }
}
