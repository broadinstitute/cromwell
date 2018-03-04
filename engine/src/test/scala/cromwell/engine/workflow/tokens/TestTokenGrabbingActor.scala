package cromwell.engine.workflow.tokens

import akka.actor.{Actor, ActorRef, Props, SupervisorStrategy}
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor.{JobExecutionTokenDispensed, JobExecutionTokenRequest}
import cromwell.util.AkkaTestUtil

/**
  * Grabs a token and doesn't let it go!
  */
class TestTokenGrabbingActor(tokenDispenser: ActorRef, tokenType: JobExecutionTokenType) extends Actor {

  var hasToken: Boolean = false

  def receive = {
    case JobExecutionTokenDispensed => hasToken = true
    case AkkaTestUtil.ThrowException => throw new RuntimeException("Test exception (don't be scared by the stack trace, it's deliberate!)")
    case AkkaTestUtil.InternalStop => context.stop(self)
  }

  tokenDispenser ! JobExecutionTokenRequest(tokenType)
}

object TestTokenGrabbingActor {

  def props(tokenDispenserActor: ActorRef, tokenType: JobExecutionTokenType) = Props(new TestTokenGrabbingActor(tokenDispenserActor, tokenType))

  class StoppingSupervisor extends Actor {
    override val supervisorStrategy = SupervisorStrategy.stoppingStrategy
    override def receive = Actor.emptyBehavior
  }
}
