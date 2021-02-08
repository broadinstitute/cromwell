package cromwell.backend.google.pipelines.common.api.clients

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.backend.google.pipelines.common.PapiInstrumentation
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.{PAPIAbortRequest, PipelinesApiAbortQueryFailed}
import cromwell.backend.google.pipelines.common.api.clients.PipelinesApiAbortClient.{PAPIAbortRequestSuccessful, PAPIOperationAlreadyCancelled, PAPIOperationHasAlreadyFinished}
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.WorkflowId
import cromwell.core.logging.JobLogging

object PipelinesApiAbortClient {
  sealed trait PAPIAbortRequestSuccess
  case class PAPIAbortRequestSuccessful(operationId: String) extends PAPIAbortRequestSuccess
  case class PAPIOperationAlreadyCancelled(operationId: String) extends PAPIAbortRequestSuccess
  case class PAPIOperationHasAlreadyFinished(operationId: String) extends PAPIAbortRequestSuccess
}

trait PipelinesApiAbortClient { this: Actor with ActorLogging with JobLogging with PapiInstrumentation =>
  def workflowId: WorkflowId

  val papiApiActor: ActorRef
  val requestFactory: PipelinesApiRequestFactory
  
  def abortJob(jobId: StandardAsyncJob) = {
    papiApiActor ! PAPIAbortRequest(workflowId, self, requestFactory.cancelRequest(jobId), jobId)
  }

  def abortActorClientReceive: Actor.Receive = {
    /* When we deprecate PAPIv1: v2+ do not seem to distinguish between cancelled and finished, unify into "terminated" or similar [PROD-444] */

    case PAPIAbortRequestSuccessful(jobId) =>
      abortSuccess()
      jobLogger.info(s"Successfully requested cancellation of $jobId")
    // In this case we could immediately return an aborted handle and spare ourselves a round of polling
    case PAPIOperationAlreadyCancelled(jobId) =>
      jobLogger.info(s"Operation $jobId was already cancelled")
    // In this case we could immediately return an aborted handle and spare ourselves a round of polling
    case PAPIOperationHasAlreadyFinished(jobId) =>
      jobLogger.info(s"Operation $jobId has already finished")
    case PipelinesApiAbortQueryFailed(jobId, e) =>
      jobLogger.error(s"Could not request cancellation of job $jobId", e)
  }
}
