package cromwell.engine.workflow.tokens

import akka.actor.{Actor, ActorRef, Props, SupervisorStrategy}
import cromwell.core.HogGroup
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor.{JobExecutionTokenDispensed, JobExecutionTokenRequest}
import cromwell.util.AkkaTestUtil.DeathTestActor

/**
  * Grabs a token and doesn't let it go!
  */
class TestTokenGrabbingActor(tokenDispenser: ActorRef, tokenType: JobExecutionTokenType) extends DeathTestActor {

  var hasToken: Boolean = false

  override def receive = stoppingReceive orElse {
    case JobExecutionTokenDispensed => hasToken = true
  }

  tokenDispenser ! JobExecutionTokenRequest(HogGroup("hogGroupA"), tokenType)
}

object TestTokenGrabbingActor {

  def props(tokenDispenserActor: ActorRef, tokenType: JobExecutionTokenType) = Props(new TestTokenGrabbingActor(tokenDispenserActor, tokenType))

  class StoppingSupervisor extends Actor {
    override val supervisorStrategy = SupervisorStrategy.stoppingStrategy
    override def receive = Actor.emptyBehavior
  }
}
