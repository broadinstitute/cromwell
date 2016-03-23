package cromwell.backend

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.event.LoggingReceive
import cromwell.backend.BackendActor._
import cromwell.backend.model.{JobDescriptor, WorkflowDescriptor}

object BackendActor {

  // Commands
  sealed trait BackendActorMessage
  case object BeforeAll extends BackendActorMessage
  case class CreateJobExecutorActor(jobDescriptor: JobDescriptor) extends BackendActorMessage
  case object AfterAll extends BackendActorMessage
  case object Validate extends BackendActorMessage

  // Events
  sealed trait BackendActorEvent

  sealed trait BeforeAllEvent extends BackendActorEvent
  case object BeforeAllSucceeded extends BeforeAllEvent
  case class BeforeAllFailed(throwable: Throwable) extends BeforeAllEvent

  sealed trait JobExecutorCreationEvent extends BackendActorEvent
  case class JobExecutorCreationSucceeded(backendJobExecutorActor: ActorRef) extends JobExecutorCreationEvent
  case class JobExecutorCreationFailed(throwable: Throwable) extends JobExecutorCreationEvent

  sealed trait AfterAllEvent extends BackendActorEvent
  case object AfterAllSucceeded extends AfterAllEvent
  case class AfterAllFailed(throwable: Throwable) extends AfterAllEvent

  sealed trait ValidationEvent extends BackendActorEvent
  case object ValidationSucceeded extends ValidationEvent
  case class ValidationFailed(throwable: Throwable) extends ValidationEvent

}

/**
  * Defines basic structure and functionality to initialize and make use of a backend through an Akka actor system.
  * Backend functions should be implemented by each custom backend.
  */
trait BackendActor extends Actor with ActorLogging {
  /**
    * Defines needed data to be able to execute a workflow.
    */
  val workflowDescriptor: WorkflowDescriptor

  def receive: Receive = LoggingReceive {
    case BeforeAll =>
      val sndr = sender()
      sndr ! beforeAll()
    case CreateJobExecutorActor(jobDescriptor: JobDescriptor) =>
      val sndr = sender()
      sndr ! createJobExecutorActor(jobDescriptor)
    case AfterAll =>
      val sndr = sender()
      sndr ! afterAll()
    case Validate => //TODO: to be discussed.
      val sndr = sender()
      sndr ! validate()
  }

  /**
    * Registers code to be executed before the backend is ready for executing jobs for the specific workflow.
    */
  def beforeAll(): BeforeAllEvent

  /**
    * Factory for creating Jobs in the specific backend.
    *
    * @param jobDescriptor All information needed to execute a task in the backend.
    * @return A BackendJobExecutor with functionality to handle the life cycle management of the task.
    */
  def createJobExecutorActor(jobDescriptor: JobDescriptor): JobExecutorCreationEvent

  /**
    * Registers code to be executed after the backend finished executing all related tasks for the specific workflow.
    */
  def afterAll(): AfterAllEvent

  /**
    * Executes validation on workflow descriptor in order to see if the workflow can be executed by the backend.
    */
  def validate(): ValidationEvent
}
