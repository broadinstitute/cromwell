package cromwell.engine

import akka.pattern.pipe
import akka.actor.{ActorRef, LoggingFSM, PoisonPill}
import cromwell.logging.WorkflowLogger
import scala.concurrent.ExecutionContext.Implicits.global

import scala.concurrent.Future

trait CromwellFSM[S, D] extends LoggingFSM[S, D] with CromwellActor {
  case class ActorFailure(failureContext: String, failure: Throwable)

  def logger: WorkflowLogger
  def sendFailureTo: Option[ActorRef]
  def failureState: S

  def fail(asyncFailure: ActorFailure, message: Any) = {
      logger.error(asyncFailure.failureContext, asyncFailure.failure)
      sendFailureTo foreach { _ ! message }
      self ! PoisonPill
      goto(failureState)
  }

  def chain(future: Future[_], nextState: S, failureContext: String) = {
    val result = future recover {
      case e => ActorFailure(failureContext, e)
    }
    pipe(result) to self

    goto(nextState)
  }
}
