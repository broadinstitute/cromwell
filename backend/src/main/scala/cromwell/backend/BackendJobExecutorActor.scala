package cromwell.backend

import akka.actor.{Actor, ActorLogging}
import akka.event.LoggingReceive
import cromwell.backend.BackendJobExecutorActor.{CleanUp, Execute, Prepare, Stop}

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
  case object PrepareCompleted extends BackendJobExecutorActorEvent
  case class PrepareFailed(throwable: Throwable) extends BackendJobExecutorActorEvent
  case object ExecutionStarted extends BackendJobExecutorActorEvent
  case class ExecutionCompleted(result: String) extends BackendJobExecutorActorEvent
  case class ExecutionFailed(throwable: Throwable) extends BackendJobExecutorActorEvent
  case object StopCompleted extends BackendJobExecutorActorEvent
  case class StopFailed(throwable: Throwable) extends BackendJobExecutorActorEvent
  case object CleanUpCompleted extends BackendJobExecutorActorEvent
  case class CleanUpFailed(throwable: Throwable) extends BackendJobExecutorActorEvent

}

/**
  * Defines basic structure and functionality to execute a job in a backend through an Akka actor system.
  * BackendJobExecutor functions should be implemented by each custom backend.
  */
trait BackendJobExecutorActor extends BackendJobExecutor with Actor with ActorLogging {
  def receive: Receive = LoggingReceive {
    case Prepare => prepare()
    case Execute => execute()
    case Stop => stop()
    case CleanUp => cleanUp()
  }
}
