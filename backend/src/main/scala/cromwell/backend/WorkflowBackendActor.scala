package cromwell.backend

import akka.actor.{ActorLogging, Actor}
import akka.event.LoggingReceive
import cromwell.backend.WorkflowBackendActor._
import cromwell.core.CallOutputs

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}
import WorkflowBackendActor.ExecutionFailed

object WorkflowBackendActor {

  // Commands
  sealed trait WorkflowBackendActorMessage

  case class Execute(jobDescriptor: BackendJobDescriptor) extends WorkflowBackendActorMessage
  case class Recover(jobDescriptor: BackendJobDescriptor) extends WorkflowBackendActorMessage
  case class Abort(jobKey: BackendJobDescriptorKey) extends WorkflowBackendActorMessage

  // Responses
  sealed trait WorkflowBackendActorResponse

  sealed trait ExecutionResponse extends WorkflowBackendActorResponse
  case class ExecutionSucceeded(jobKey: BackendJobDescriptorKey, callOutputs: CallOutputs) extends ExecutionResponse
  case class ExecutionFailed(jobKey: BackendJobDescriptorKey, throwable: Throwable) extends ExecutionResponse
  case class ExecutionFailedRetryable(jobKey: BackendJobDescriptorKey, throwable: Throwable) extends ExecutionResponse

  sealed trait AbortResponse extends WorkflowBackendActorResponse
  case class AbortSucceeded(jobKey: BackendJobDescriptorKey) extends AbortResponse
  case class AbortFailed(jobKey: BackendJobDescriptorKey, throwable: Throwable) extends AbortResponse
}

/**
  * Workflow-level actor for executing, recovering and aborting jobs.
  */
trait WorkflowBackendActor extends Actor with ActorLogging {

  protected implicit def ec: ExecutionContext
  protected def workflowDescriptor: BackendWorkflowDescriptor
  protected def configurationDescriptor: BackendConfigurationDescriptor

  def receive: Receive = LoggingReceive {
    case Execute(jobDescriptor) => performActionThenRespond(execute, jobDescriptor, onFailure = descriptorToExecutionFailure)
    case Recover(jobDescriptor) => performActionThenRespond(recover, jobDescriptor, onFailure = descriptorToExecutionFailure)
    case Abort(jobKey)          => performActionThenRespond(abort, jobKey, onFailure = AbortFailed)
  }

  private def performActionThenRespond[A](operation: A => Future[WorkflowBackendActorResponse],
                                          operationArgument: A,
                                          onFailure: (A, Throwable) => WorkflowBackendActorResponse) = {
    // To avoid sending a response to the wrong `sender` if more than one `execute` arrives before the first one finishes,
    // for example:
    val sndr = sender()
    operation(operationArgument) onComplete {
      case Success(r) => sndr ! r
      case Failure(t) => sndr ! onFailure(operationArgument, t)
    }
  }

  // We need this for receive because we can't do `onFailure = ExecutionFailure` directly - because BackendJobDescriptor =/= BackendJobDescriptorKey
  private def descriptorToExecutionFailure = (jd: BackendJobDescriptor, t: Throwable) => ExecutionFailed(jd.key, t)

  /**
    * Execute a new job.
    */
  def execute(jobDescriptor: BackendJobDescriptor): Future[ExecutionResponse]

  /**
    * Restart or resume a previously-started job.
    */
  def recover(jobDescriptor: BackendJobDescriptor): Future[ExecutionResponse]

  /**
    * Abort a running job.
    */
  def abort(jobKey: BackendJobDescriptorKey): Future[AbortResponse]
}
