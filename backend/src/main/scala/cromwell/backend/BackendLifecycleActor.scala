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
  case object AbortWorkflowCommand extends BackendWorkflowLifecycleActorCommand
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
                                         onFailure: (Throwable) => BackendWorkflowLifecycleActorResponse) = {
    val respondTo: ActorRef = sender
    operation onComplete {
      case Success(r) => respondTo ! r
      case Failure(t) => respondTo ! onFailure(t)
    }
  }
}

trait BackendWorkflowLifecycleActor extends BackendLifecycleActor {
  /**
    * The workflow descriptor for the workflow in which this Backend is being used
    */
  protected def workflowDescriptor: BackendWorkflowDescriptor

  /**
    * The subset of calls which this backend will be expected to run
    */
  protected def calls: Seq[Call]
}

trait BackendJobLifecycleActor extends BackendLifecycleActor {
  protected def jobDescriptor: BackendJobDescriptor
}
