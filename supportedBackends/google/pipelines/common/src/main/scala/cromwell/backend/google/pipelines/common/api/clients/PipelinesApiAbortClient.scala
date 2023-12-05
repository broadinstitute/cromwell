package cromwell.backend.google.pipelines.common.api.clients

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.backend.google.pipelines.common.PapiInstrumentation
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.{
  PAPIAbortRequest,
  PipelinesApiAbortQueryFailed
}
import cromwell.backend.google.pipelines.common.api.clients.PipelinesApiAbortClient.{
  PAPIAbortRequestSuccessful,
  PAPIOperationIsAlreadyTerminal
}
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.WorkflowId
import cromwell.core.logging.JobLogging

object PipelinesApiAbortClient {
  sealed trait PAPIAbortRequestSuccess
  case class PAPIAbortRequestSuccessful(operationId: String) extends PAPIAbortRequestSuccess
  // The operation is no longer running. Maybe it was already cancelled, maybe it finished on its own. We don't know
  // the details and for abort they don't really matter.
  case class PAPIOperationIsAlreadyTerminal(operationId: String) extends PAPIAbortRequestSuccess
}

trait PipelinesApiAbortClient { this: Actor with ActorLogging with JobLogging with PapiInstrumentation =>
  def workflowId: WorkflowId

  val papiApiActor: ActorRef
  val requestFactory: PipelinesApiRequestFactory

  def abortJob(jobId: StandardAsyncJob) =
    papiApiActor ! PAPIAbortRequest(workflowId, self, requestFactory.cancelRequest(jobId), jobId)

  def abortActorClientReceive: Actor.Receive = {
    case PAPIAbortRequestSuccessful(jobId) =>
      abortSuccess()
      jobLogger.info(s"Successfully requested cancellation of $jobId")
    // In this case we could immediately return an aborted handle and spare ourselves a round of polling
    case PAPIOperationIsAlreadyTerminal(jobId) =>
      jobLogger.info(s"Operation $jobId has already finished")
    case PipelinesApiAbortQueryFailed(jobId, e) =>
      jobLogger.error(s"Could not request cancellation of job $jobId", e)
  }
}
