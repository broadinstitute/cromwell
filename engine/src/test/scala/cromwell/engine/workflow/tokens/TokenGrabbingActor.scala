package cromwell.engine.workflow.tokens

import akka.actor.{Actor, ActorRef, Props, SupervisorStrategy}
import cromwell.core.JobExecutionToken
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor.{JobExecutionTokenDenied, JobExecutionTokenDispensed, JobExecutionTokenRequest}
import cromwell.engine.workflow.tokens.TokenGrabbingActor.{InternalStop, ThrowException}

/**
  * Grabs a token and doesn't let it go!
  */
class TokenGrabbingActor(tokenDispenser: ActorRef, tokenType: JobExecutionTokenType) extends Actor {

  var token: Option[JobExecutionToken] = None
  var rejections = 0

  def receive = {
    case JobExecutionTokenDispensed(dispensedToken) => token = Option(dispensedToken)
    case JobExecutionTokenDenied(positionInQueue) => rejections += 1
    case ThrowException => throw new RuntimeException("Test exception (don't be scared by the stack trace, it's deliberate!)")
    case InternalStop => context.stop(self)
  }

  tokenDispenser ! JobExecutionTokenRequest(tokenType)
}

object TokenGrabbingActor {

  def props(tokenDispenserActor: ActorRef, tokenType: JobExecutionTokenType) = Props(new TokenGrabbingActor(tokenDispenserActor, tokenType))

  case object ThrowException
  case object InternalStop

  class StoppingSupervisor extends Actor {
    override val supervisorStrategy = SupervisorStrategy.stoppingStrategy
    override def receive = Actor.emptyBehavior
  }
}
