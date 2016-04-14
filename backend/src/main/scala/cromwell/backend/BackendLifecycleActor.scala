package cromwell.backend

import akka.actor.{ActorRef, Actor}
import cromwell.backend.BackendLifecycleActor._
import wdl4s.Call

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

object BackendLifecycleActor {

  trait BackendWorkflowLifecycleActorMessage

  /*
   * Commands
   */
  trait BackendWorkflowLifecycleActorCommand extends BackendWorkflowLifecycleActorMessage
  case object AbortWorkflow extends BackendWorkflowLifecycleActorCommand
  final case class AbortJob(jobKey: BackendJobDescriptorKey) extends BackendWorkflowLifecycleActorCommand

  /*
   * Responses
   */
  trait BackendWorkflowLifecycleActorResponse extends BackendWorkflowLifecycleActorMessage

  sealed trait JobAbortResponse extends BackendWorkflowLifecycleActorResponse
  sealed trait WorkflowAbortResponse extends BackendWorkflowLifecycleActorResponse

  case object BackendWorkflowAbortSucceededResponse extends WorkflowAbortResponse
  final case class BackendWorkflowAbortFailedResponse(throwable: Throwable) extends WorkflowAbortResponse
  final case class BackendJobExecutionAbortSucceededResponse(jobKey: BackendJobDescriptorKey) extends JobAbortResponse
  final case class BackendJobExecutionAbortFailedResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable) extends JobAbortResponse
}

trait BackendLifecycleActor extends Actor {

  /**
    * The execution context for the actor
    */
  protected implicit def ec: ExecutionContext = context.dispatcher

  /**
    * The set of calls which this backend will be expected to run
    */
  protected def calls: Seq[Call]

  /**
    * The workflow descriptor for the workflow in which this Backend is being used
    */
  protected def workflowDescriptor: BackendWorkflowDescriptor

  /**
    * The configuration for the backend, in the context of the entire Cromwell configuration file.
    */
  protected def configurationDescriptor: BackendConfigurationDescriptor

  protected def performActionThenRespond(operation: => Future[BackendWorkflowLifecycleActorResponse],
                                         onFailure: (Throwable) => BackendWorkflowLifecycleActorResponse) = {
    val respondTo: ActorRef = sender
    operation onComplete {
      case Success(r) => respondTo ! r
      case Failure(t) => respondTo ! onFailure(t)
    }
  }
}
