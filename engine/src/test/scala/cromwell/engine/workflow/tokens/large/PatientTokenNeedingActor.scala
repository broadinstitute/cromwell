package cromwell.engine.workflow.tokens.large

import akka.actor.{Actor, ActorRef, Props}
import cromwell.core.HogGroup
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor.{JobExecutionTokenDispensed, JobExecutionTokenRequest, JobExecutionTokenReturn}
import cromwell.engine.workflow.tokens.large.PatientTokenNeedingActor.{AllDone, Begin, ImBusy, RequestToken}
import org.joda.time.DateTime

import scala.concurrent.duration._
import scala.util.Random

/**
  * The sort-of equivalent of an EngineJobExecutionActor -
  *
  * I need a token in order to start working, and don't know anything about how many others like me there are also
  * wanting tokens.
  *
  * I'm happy to wait as long as necessary to get my token.
  */
class PatientTokenNeedingActor(tokenDispenser: ActorRef, tokenType: JobExecutionTokenType, hogGroup: String) extends Actor {

  var hasToken: Boolean = false
  var starter: ActorRef = _
  var startTime: DateTime = _

  /**
    * When told to begin, I schedule myself to ask for a token after a random pause
    * When I get a token, I "work" by scheduling an AllDone to be sent back to myself after a longer pause
    * When I get the AllDone, I forward that status on to whoever started me, then return the token
    */
  override def receive = {
    case Begin =>
      starter = sender
      val requestDelay = Random.nextInt(30) + 1
      context.system.scheduler.scheduleOnce(requestDelay.millis, self, RequestToken)(context.dispatcher)
      ()
    case RequestToken =>
      tokenDispenser ! JobExecutionTokenRequest(HogGroup(hogGroup), tokenType)
      startTime = DateTime.now()
    case JobExecutionTokenDispensed =>
      context.system.scheduler.scheduleOnce(1.seconds, self, AllDone)(context.dispatcher)
      starter ! ImBusy(DateTime.now().getMillis - startTime.getMillis)
    case AllDone =>
      starter ! AllDone
      tokenDispenser ! JobExecutionTokenReturn
      context.stop(self)
    case other =>
      throw new Exception(s"Bad message received: $other")
  }
}

object PatientTokenNeedingActor {
  // Tells this actor to begin
  case object Begin
  // Tells this actor to request a token
  private case object RequestToken
  // Indicate to parent that I'm busy
  case class ImBusy(msInQueue: Long)
  // Indicate to myself that I'm done (and gets forwarded to my parent)
  case object AllDone

  def props(tokenDispenser: ActorRef, tokenType: JobExecutionTokenType, hogGroup: String): Props = Props(new PatientTokenNeedingActor(tokenDispenser, tokenType, hogGroup))
}
