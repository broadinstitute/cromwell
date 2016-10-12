package cromwell.backend

import akka.actor.{Actor, ActorRef}
import cromwell.backend.BackendLifecycleActor._
import cromwell.core.logging.{JobLogging, WorkflowLogging}
import wdl4s.TaskCall

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

object BackendLifecycleActor {

  trait BackendWorkflowLifecycleActorMessage

  /*
   * Commands
   */
  trait BackendWorkflowLifecycleActorCommand extends BackendWorkflowLifecycleActorMessage
  case object AbortJobCommand extends BackendWorkflowLifecycleActorCommand

  /*
   * Responses
   */
  trait BackendWorkflowLifecycleActorResponse extends BackendWorkflowLifecycleActorMessage
  case object BackendActorAbortedResponse extends BackendWorkflowLifecycleActorResponse
}

trait BackendLifecycleActor extends Actor {

  /**
    * The execution context for the actor
    */
  protected implicit def ec: ExecutionContext = context.dispatcher

  /**
    * The configuration for the backend, in the context of the entire Cromwell configuration file.
    */
  protected def configurationDescriptor: BackendConfigurationDescriptor

  protected def performActionThenRespond(operation: => Future[BackendWorkflowLifecycleActorResponse],
                                         onFailure: (Throwable) => BackendWorkflowLifecycleActorResponse,
                                         andThen: => Unit = ()) = {
    val respondTo: ActorRef = sender
    operation onComplete {
      case Success(r) =>
        respondTo ! r
        andThen
      case Failure(t) =>
        respondTo ! onFailure(t)
        andThen
    }
  }
}

trait BackendWorkflowLifecycleActor extends BackendLifecycleActor with WorkflowLogging {

  //For Logging and boilerplate
  override lazy final val workflowIdForLogging = workflowDescriptor.id

  /**
    * The workflow descriptor for the workflow in which this Backend is being used
    */
  protected def workflowDescriptor: BackendWorkflowDescriptor

  /**
    * The subset of calls which this backend will be expected to run
    */
  protected def calls: Set[TaskCall]
}

trait BackendJobLifecycleActor extends BackendLifecycleActor with JobLogging {
  //For Logging and boilerplate
  override lazy val workflowId = jobDescriptor.workflowDescriptor.id
  override lazy val jobTag = jobDescriptor.key.tag

  protected def jobDescriptor: BackendJobDescriptor
}
