package cromwell.engine.workflow.tokens

import java.util.UUID

import cromwell.core.JobExecutionToken
import JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.TokenPool.TokenPoolPop

import scala.language.postfixOps

trait TokenPool {
  def currentLoans: Set[JobExecutionToken]
  def pop(): TokenPoolPop
  def push(jobExecutionToken: JobExecutionToken): TokenPool
}

object TokenPool {

  case class TokenPoolPop(newTokenPool: TokenPool, poppedItem: Option[JobExecutionToken])

  def apply(tokenType: JobExecutionTokenType): TokenPool = {
    tokenType.maxPoolSize map { ps =>
      val pool = (1 to ps toList) map { _ => JobExecutionToken(tokenType, UUID.randomUUID()) }
      SizedTokenPool(pool, Set.empty)
    } getOrElse {
      InfiniteTokenPool(tokenType, Set.empty)
    }
  }
}

final case class SizedTokenPool(pool: List[JobExecutionToken], override val currentLoans: Set[JobExecutionToken]) extends TokenPool {

  override def pop(): TokenPoolPop = pool match {
    case head :: tail => TokenPoolPop(SizedTokenPool(tail, currentLoans + head), Option(head))
    case Nil => TokenPoolPop(SizedTokenPool(List.empty, currentLoans), None)
  }


  override def push(token: JobExecutionToken): SizedTokenPool = {
    if (currentLoans.contains(token)) { SizedTokenPool(pool :+ token, currentLoans - token) } else this
  }
}

final case class InfiniteTokenPool(tokenType: JobExecutionTokenType, override val currentLoans: Set[JobExecutionToken]) extends TokenPool {
  override def pop() = {
    val newToken = JobExecutionToken(tokenType, UUID.randomUUID())
    TokenPoolPop(InfiniteTokenPool(tokenType, currentLoans + newToken), Option(newToken))
  }
  override def push(token: JobExecutionToken): InfiniteTokenPool = if (currentLoans.contains(token)) { InfiniteTokenPool(tokenType, currentLoans - token) } else this
}
