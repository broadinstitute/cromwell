package cromwell.backend

import akka.actor.{Actor, ActorLogging}
import akka.event.LoggingReceive
import cromwell.backend.BackendJobExecutorActor._
import cromwell.backend.model.JobDescriptor

object BackendJobExecutorActor {

  // Commands
  sealed trait BackendJobExecutorActorMessage
  case object Prepare extends BackendJobExecutorActorMessage
  case object Execute extends BackendJobExecutorActorMessage
  case object Stop extends BackendJobExecutorActorMessage
  case object CleanUp extends BackendJobExecutorActorMessage

  //TODO: identify other generic events within job execution functionality.
  // Events
  sealed trait BackendJobExecutorActorEvent

  sealed trait PrepareEvent extends BackendJobExecutorActorEvent
  case object PrepareSucceeded extends PrepareEvent
  case class PrepareFailed(throwable: Throwable) extends PrepareEvent

  sealed trait ExecutionEvent extends BackendJobExecutorActorEvent
  case class ExecutionSucceeded(result: String) extends ExecutionEvent
  case class ExecutionFailed(throwable: Throwable) extends ExecutionEvent

  sealed trait StopEvent extends BackendJobExecutorActorEvent
  case object StopSucceeded extends StopEvent
  case class StopFailed(throwable: Throwable) extends StopEvent

  sealed trait CleanUpEvent extends BackendJobExecutorActorEvent
  case object CleanUpSucceeded extends CleanUpEvent
  case class CleanUpFailed(throwable: Throwable) extends CleanUpEvent

}

/**
  * Defines basic structure and functionality to execute a job in a backend through an Akka actor system.
  * BackendJobExecutor functions should be implemented by each custom backend.
  */
trait BackendJobExecutorActor extends Actor with ActorLogging {
  /**
    * Defines needed data to be able to execute a job.
    */
  val jobDescriptor: JobDescriptor

  def receive: Receive = LoggingReceive {
    case Prepare =>
      val sndr = sender()
      sndr ! prepare()
    case Execute =>
      val sndr = sender()
      sndr ! execute()
    case Stop =>
      val sndr = sender()
      sndr ! stop()
    case CleanUp =>
      val sndr = sender()
      sndr ! cleanUp()
  }

  /**
    * Prepare the task and context for execution.
    */
  def prepare(): PrepareEvent

  /**
    * Executes task in given context.
    */
  def execute(): ExecutionEvent

  /**
    * Stops a task execution.
    */
  def stop(): StopEvent

  /**
    * Performs a cleanUp after the task was executed.
    */
  def cleanUp(): CleanUpEvent
}
