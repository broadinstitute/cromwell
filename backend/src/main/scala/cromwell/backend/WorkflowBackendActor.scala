package cromwell.backend

import akka.actor.{ActorLogging, Actor}
import akka.event.LoggingReceive
import cromwell.backend.WorkflowBackendActor._
import cromwell.core.CallOutputs
import wdl4s.Call

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

  sealed trait ValidationResponse extends WorkflowBackendActorResponse
  case object ValidationSuccess extends ValidationResponse
  case class ValidationFailed(reason: Throwable) extends ValidationResponse
}

/**
  * Workflow-level actor for executing, recovering and aborting jobs.
  */
trait WorkflowBackendActor extends Actor with ActorLogging {

  protected implicit def ec: ExecutionContext
  protected def calls: Seq[Call]
  protected def workflowDescriptor: BackendWorkflowDescriptor
  protected def configurationDescriptor: BackendConfigurationDescriptor

  def receive: Receive = LoggingReceive {
    case Execute(jobDescriptor) => performActionThenRespond( () => execute(jobDescriptor), onFailure = executionFailed(jobDescriptor.key))
    case Recover(jobDescriptor) => performActionThenRespond( () => recover(jobDescriptor), onFailure = executionFailed(jobDescriptor.key))
    case Abort(jobKey)          => performActionThenRespond( () => abort(jobKey), onFailure = abortFailed(jobKey))
  }

  override def preStart = {
    // Do everything necessary before this actor is ready. Currently just validation:
    performActionThenRespond(validate, onFailure = ValidationFailed)
  }

  private def performActionThenRespond(operation: () => Future[WorkflowBackendActorResponse],
                                       onFailure: (Throwable) => WorkflowBackendActorResponse) = {
    // To avoid sending a response to the wrong `sender` if more than one `execute` arrives before the first one finishes,
    // for example:
    val sndr = sender()
    operation.apply onComplete {
      case Success(r) => sndr ! r
      case Failure(t) => sndr ! onFailure(t)
    }
  }

  // We need this for receive because we can't do `onFailure = ExecutionFailure` directly - because BackendJobDescriptor =/= BackendJobDescriptorKey
  private def executionFailed(key: BackendJobDescriptorKey) = (t: Throwable) => ExecutionFailed(key, t)
  private def abortFailed(key: BackendJobDescriptorKey) = (t: Throwable) => AbortFailed(key, t)

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

  /**
    * Validate that this WorkflowBackendActor can run all of the calls that it's been assigned
    */
  def validate(): Future[ValidationResponse]
}
