package cromwell.backend

import akka.actor.ActorLogging
import akka.event.LoggingReceive
import cromwell.backend.BackendLifecycleActor._
import cromwell.backend.BackendWorkflowInitializationActor._

import scala.concurrent.Future

object BackendWorkflowInitializationActor {

  // Commands
  sealed trait BackendWorkflowInitializationActorCommand extends BackendWorkflowLifecycleActorCommand
  case object Initialize extends BackendWorkflowInitializationActorCommand
  case object Abort extends BackendWorkflowInitializationActorCommand

  // Responses
  sealed trait BackendWorkflowInitializationActorResponse extends BackendWorkflowLifecycleActorResponse
  sealed trait InitializationResponse extends BackendWorkflowInitializationActorResponse
  case object InitializationSuccess extends InitializationResponse
  case class InitializationFailed(reason: Throwable) extends Exception with InitializationResponse

}

/**
  * Workflow-level actor for executing, recovering and aborting jobs.
  */
trait BackendWorkflowInitializationActor extends BackendWorkflowLifecycleActor with ActorLogging {

  def receive: Receive = LoggingReceive {
    case Initialize => performActionThenRespond(initSequence(), onFailure = InitializationFailed)
    case AbortWorkflowCommand => abortInitialization()
  }

  /**
    * Our predefined sequence to run during preStart
    */
  final def initSequence() = for {
    _ <- validate()
    _ <- beforeAll()
  } yield InitializationSuccess

  /**
    * Abort all initializations.
    */
  def abortInitialization(): Unit

  /**
    * A call which happens before anything else runs
    */
  def beforeAll(): Future[Unit]

  /**
    * Validate that this WorkflowBackendActor can run all of the calls that it's been assigned
    */
  def validate(): Future[Unit]

}
