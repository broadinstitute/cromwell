package cromwell.backend

import akka.actor.ActorLogging
import akka.event.LoggingReceive
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.BackendLifecycleActor._
import cromwell.core.CallOutputs

import scala.concurrent.Future
import BackendJobExecutionActor.BackendJobExecutionFailedResponse

object BackendJobExecutionActor {

  // Commands
  sealed trait BackendJobExecutionActorCommand extends BackendWorkflowLifecycleActorCommand
  case object ExecuteJobCommand extends BackendJobExecutionActorCommand
  case object RecoverJobCommand extends BackendJobExecutionActorCommand

  // Responses
  sealed trait BackendJobExecutionActorResponse extends BackendWorkflowLifecycleActorResponse

  sealed trait BackendJobExecutionResponse extends BackendJobExecutionActorResponse
  case class BackendJobExecutionSucceededResponse(jobKey: BackendJobDescriptorKey, callOutputs: CallOutputs) extends BackendJobExecutionResponse
  case class BackendJobExecutionFailedResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable) extends BackendJobExecutionResponse
  case class BackendJobExecutionFailedRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable) extends BackendJobExecutionResponse
}

/**
  * Workflow-level actor for executing, recovering and aborting jobs.
  */
trait BackendJobExecutionActor extends BackendJobLifecycleActor with ActorLogging {

  def receive: Receive = LoggingReceive {
    case ExecuteJobCommand => performActionThenRespond(execute, onFailure = executionFailed)
    case RecoverJobCommand => performActionThenRespond(recover, onFailure = executionFailed)
    case AbortJob          => performActionThenRespond(abortJob, onFailure = abortFailed)
  }

  // We need this for receive because we can't do `onFailure = ExecutionFailure` directly - because BackendJobDescriptor =/= BackendJobDescriptorKey
  private def executionFailed = (t: Throwable) => BackendJobExecutionFailedResponse(jobDescriptor.key, t)
  private def abortFailed = (t: Throwable) => BackendJobExecutionAbortFailedResponse(jobDescriptor.key, t)

  /**
    * Execute a new job.
    */
  def execute: Future[BackendJobExecutionResponse]

  /**
    * Restart or resume a previously-started job.
    */
  def recover: Future[BackendJobExecutionResponse]

  /**
    * Abort a running job.
    */
  def abortJob: Future[JobAbortResponse]
}
