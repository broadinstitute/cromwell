package cromwell.engine.workflow.tokens

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.core.JobExecutionToken
import JobExecutionToken._
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor._

import scala.language.postfixOps

class JobExecutionTokenDispenserActor extends Actor with ActorLogging {

  /**
    * Lazily created token pool. We only create a pool for a token type when we need it
    */
  var tokenPools: Map[JobExecutionTokenType, TokenPool] = Map.empty

  // TODO: Record which actors are given which token, and reclaim them on terminate.

  override def receive = {
    case JobExecutionTokenRequest(tokenType) =>
      val pop = tokenPools.getOrElse(tokenType, TokenPool(tokenType)).pop()
      tokenPools += tokenType -> pop.newTokenPool
      sender ! (pop.poppedItem map { token =>
        JobExecutionTokenDispensed(token)
      } getOrElse { JobExecutionTokenDenied })
    case JobExecutionTokenReturn(token) =>
      tokenPools.get(token.jobExecutionTokenType) foreach { pool =>
        tokenPools += token.jobExecutionTokenType -> pool.push(token)
      }
  }
}

object JobExecutionTokenDispenserActor {

  def props = Props(new JobExecutionTokenDispenserActor)

  case class JobExecutionTokenRequest(jobExecutionTokenType: JobExecutionTokenType)
  case class JobExecutionTokenReturn(jobExecutionToken: JobExecutionToken)

  // TODO: Implement a queue based system!
  case object JobExecutionTokenDenied // TODO: (positionInQueue: Integer)
  case class JobExecutionTokenDispensed(jobExecutionToken: JobExecutionToken)

}