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
  final case class ExecuteJobCommand(jobDescriptor: BackendJobDescriptor) extends BackendJobExecutionActorCommand
  final case class RecoverJobCommand(jobDescriptor: BackendJobDescriptor) extends BackendJobExecutionActorCommand

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
    case ExecuteJobCommand(jobDescriptor) => performActionThenRespond(execute(jobDescriptor), onFailure = executionFailed(jobDescriptor.key))
    case RecoverJobCommand(jobDescriptor) => performActionThenRespond(recover(jobDescriptor), onFailure = executionFailed(jobDescriptor.key))
    case AbortJob(jobKey)                 => performActionThenRespond(abortJob(jobKey), onFailure = abortFailed(jobKey))
  }

  // We need this for receive because we can't do `onFailure = ExecutionFailure` directly - because BackendJobDescriptor =/= BackendJobDescriptorKey
  private def executionFailed(key: BackendJobDescriptorKey) = (t: Throwable) => BackendJobExecutionFailedResponse(key, t)
  private def abortFailed(key: BackendJobDescriptorKey) = (t: Throwable) => BackendJobExecutionAbortFailedResponse(key, t)

  override def jobDescriptor: BackendJobDescriptor

  /**
    * Execute a new job.
    */
  def execute(): Future[BackendJobExecutionResponse]

  /**
    * Restart or resume a previously-started job.
    */
  def recover(): Future[BackendJobExecutionResponse]

  /**
    * Abort a running job.
    */
  def abortJob(): Future[JobAbortResponse]
}
