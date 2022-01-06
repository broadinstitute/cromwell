package cromwell.engine.workflow.tokens

import akka.actor.{Actor, ActorRef, Props, SupervisorStrategy}
import cromwell.core.HogGroup
import cromwell.core.JobToken.JobTokenType
import cromwell.engine.workflow.tokens.JobTokenDispenserActor.{JobTokenDispensed, JobTokenRequest}
import cromwell.util.AkkaTestUtil.DeathTestActor

/**
  * Grabs a token and doesn't let it go!
  */
class TestTokenGrabbingActor(tokenDispenser: ActorRef, tokenType: JobTokenType) extends DeathTestActor {

  var hasToken: Boolean = false

  override def receive = stoppingReceive orElse {
    case JobTokenDispensed => hasToken = true
  }

  tokenDispenser ! JobTokenRequest(HogGroup("hogGroupA"), tokenType)
}

object TestTokenGrabbingActor {

  def props(tokenDispenserActor: ActorRef, tokenType: JobTokenType) = Props(new TestTokenGrabbingActor(tokenDispenserActor, tokenType))

  class StoppingSupervisor extends Actor {
    override val supervisorStrategy = SupervisorStrategy.stoppingStrategy
    override def receive = Actor.emptyBehavior
  }
}
