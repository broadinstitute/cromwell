package cromwell.engine

import akka.actor.{ActorRef, LoggingFSM, PoisonPill}
import akka.pattern.pipe
import cromwell.logging.WorkflowLogger

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future

object CromwellFSM {
  trait ActorFailure {
    def failureContext: String
    def failure: Throwable
  }
}

trait CromwellFSM[S, D] extends LoggingFSM[S, D] with CromwellActor {
  import cromwell.engine.CromwellFSM._

  def logger: WorkflowLogger
  def sendFailureTo: Option[ActorRef]
  def failureState: S

  def fail(actorFailure: ActorFailure, failureMapper: (ActorFailure => Any) = identity[ActorFailure]) = {
      logger.error(actorFailure.failureContext, actorFailure.failure)
      sendFailureTo foreach { _ ! failureMapper(actorFailure) }
      self ! PoisonPill
      goto(failureState)
  }
}

trait RetryableFSM[S, D] { this: CromwellFSM[S, D] =>
  import cromwell.engine.CromwellFSM._

  def isRetryable(failure: ActorFailure): Boolean

  def retryOrFail(failure: ActorFailure, originalMessage: Any, failureMapper: (ActorFailure => Any) = identity[ActorFailure]) = {
    if (isRetryable(failure)) {
      logger.error(failure.failureContext, failure.failure)
      logger.error("Caught retryable failure - retrying.")
      self ! originalMessage
      stay()
    } else {
      fail(failure, failureMapper)
    }
  }
}