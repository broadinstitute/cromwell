package cromwell.backend

import akka.actor.{ActorRef, Actor}
import cromwell.backend.BackendLifecycleActor._
import cromwell.core.Evaluator
import wdl4s.values.WdlValue
import wdl4s.{LocallyQualifiedName, Call}

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Try, Failure, Success}

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

  // Boilerplate code to load the configuration from the descriptor
  lazy val backendConfiguration = configurationDescriptor.config.getConfig(configurationDescriptor.configPath)

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
  protected def evaluateInputs(evaluator: Evaluator): Map[LocallyQualifiedName, Try[WdlValue]] = {
    val a = jobDescriptor.unevaluatedInputs mapValues evaluator.evaluate
    a
  }
}
