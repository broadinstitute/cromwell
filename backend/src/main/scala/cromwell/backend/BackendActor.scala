package cromwell.backend

import akka.actor.{Actor, ActorLogging}
import akka.event.LoggingReceive
import cromwell.backend.BackendActor._
import cromwell.backend.model.{JobDescriptor, WorkflowDescriptor}

import scala.concurrent.Future
import scala.util.{Success, Failure}

object BackendActor {

  // Commands
  sealed trait BackendActorMessage
  case class Execute(jobDescriptor: JobDescriptor) extends BackendActorMessage
  case class Abort(jobDescriptor: JobDescriptor) extends BackendActorMessage
  case object AbortAll extends BackendActorMessage
  case class Recover(jobDescriptor: JobDescriptor) extends BackendActorMessage

  // Events
  sealed trait BackendActorEvent

  sealed trait ExecutionEvent extends BackendActorEvent
  case class ExecutionSucceeded(result: String) extends ExecutionEvent //Result TBD.
  case class ExecutionFailed(throwable: Throwable) extends ExecutionEvent

  sealed trait AbortEvent extends BackendActorEvent
  case class AbortSucceeded(jobDescriptor: JobDescriptor) extends AbortEvent
  case class AbortFailed(throwable: Throwable) extends AbortEvent

  sealed trait AbortAllEvent extends BackendActorEvent
  case class AbortAllSucceeded(jobDescriptor: List[JobDescriptor]) extends AbortAllEvent //Attributes TBD.
  case class AbortAllFailed(throwable: Throwable) extends AbortAllEvent

  sealed trait RecoverEvent extends BackendActorEvent
  case class RecoverSucceeded(jobDescriptor: JobDescriptor) extends RecoverEvent //Same happens here.
  case class RecoverFailed(throwable: Throwable) extends RecoverEvent

  sealed trait ValidationEvent extends BackendActorEvent
  case object ValidationSucceeded extends ValidationEvent
  case class ValidationFailed(throwable: Throwable) extends ValidationEvent

  sealed trait BeforeAllEvent extends BackendActorEvent
  case object BeforeAllSucceeded extends BeforeAllEvent
  case class BeforeAllFailed(throwable: Throwable) extends BeforeAllEvent

  sealed trait AfterAllEvent extends BackendActorEvent
  case object AfterAllSucceeded extends AfterAllEvent
  case class AfterAllFailed(throwable: Throwable) extends AfterAllEvent

}

/**
  * Defines basic structure and functionality to initialize and make use of a backend through an Akka actor system.
  * Backend functions should be implemented by each custom backend.
  */
trait BackendActor extends Actor with ActorLogging {
  import context.dispatcher

  /**
    * Defines needed data to be able to execute a workflow.
    */
  val workflowDescriptor: WorkflowDescriptor

  def receive: Receive = LoggingReceive {
    case Execute(jobDescriptor: JobDescriptor) =>
      val sndr = sender()
      execute(jobDescriptor) onComplete {
        case Success(executionEvent: ExecutionEvent) =>
          sndr ! executionEvent
        case Failure(exception: Throwable) =>
          sndr ! ExecutionFailed(exception)
      }
    case Abort(jobDescriptor: JobDescriptor) =>
      sender() ! abort(jobDescriptor)
    case AbortAll =>
      sender() ! abortAll()
    case Recover =>
      sender() ! recover()
  }

  /**
    * Makes mandatory the execution of validation and beforeAll.
    * Validation: intention here is to check if runtime attributes complies with backend specifications.
    * BeforeAll: Execute any required pre-condition before making backend available for executions.
    */
  override def preStart(): Unit = {
    validate()
    beforeAll()
  }

  /**
    * Makes mandatory the execution of afterAll.
    * AfterAll: Execute any required post-condition after shutting down the backend actor.
    */
  override def postStop(): Unit = {
    afterAll()
  }

  /**
    * Executes a job base on a job description.
    *
    * @param jobDescriptor All information needed to execute a task in the backend.
    * @return An ExecutionEvent with the result.
    */
  def execute(jobDescriptor: JobDescriptor): Future[ExecutionEvent]

  /**
    * Stops a job execution.
    *
    * @param jobDescriptor All information related to a task.
    * @return An AbortEvent with the result.
    */
  def abort(jobDescriptor: JobDescriptor): AbortEvent

  /**
    * Stops all job execution.
    *
    * @return An AbortAllEvent with the result.
    */
  def abortAll(): AbortAllEvent

  /**
    * Tries to recover last status and information on an on going job in the backend.
    *
    * @return An RecoverEvent with the result.
    */
  //TODO: Need more details on this in order to define it.
  def recover(): RecoverEvent

  /**
    * Executes validation on workflow descriptor in order to see if the workflow can be executed by the backend.
    */
  def validate(): Unit

  /**
    * Registers code to be executed before the backend is ready for executing jobs for the specific workflow.
    */
  def beforeAll(): Unit

  /**
    * Registers code to be executed after the backend finished executing all related tasks for the specific workflow.
    */
  def afterAll(): Unit

}
