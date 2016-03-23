package cromwell.backend

import akka.actor.{ActorRef, Actor, ActorLogging}
import akka.event.LoggingReceive
import cromwell.backend.BackendActor.{Validate, AfterAll, BeforeAll, CreateJobExecutor}
import cromwell.backend.model.JobDescriptor

object BackendActor {

  // Commands
  sealed trait BackendActorMessage
  case object BeforeAll extends BackendActorMessage
  case class CreateJobExecutor(jobDescriptor: JobDescriptor) extends BackendActorMessage
  case object AfterAll extends BackendActorMessage
  case object Validate extends BackendActorMessage

  // Events
  sealed trait BackendActorEvent
  case object BeforeAllCompleted extends BackendActorEvent
  case class BeforeAllFailed(throwable: Throwable) extends BackendActorEvent
  case class JobExecutorCreated(backendJobExecutor: ActorRef) extends BackendActorEvent
  case class JobExecutorCreationFailed(throwable: Throwable) extends BackendActorEvent
  case object AfterAllCompleted extends BackendActorEvent
  case class AfterAllFailed(throwable: Throwable) extends BackendActorEvent
  case object ValidationCompleted extends BackendActorEvent
  case class ValidationFailed(throwable: Throwable) extends BackendActorEvent

}

/**
  * Defines basic structure and functionality to initialize and make use of a backend through an Akka actor system.
  * Backend functions should be implemented by each custom backend.
  */
trait BackendActor extends Backend with Actor with ActorLogging {
  def receive: Receive = LoggingReceive {
    case BeforeAll => beforeAll()
    case CreateJobExecutor(jobDescriptor: JobDescriptor) => createJobExecutor(jobDescriptor)
    case AfterAll => afterAll()
    case Validate => validate() //TODO: to be discussed.
  }
}
