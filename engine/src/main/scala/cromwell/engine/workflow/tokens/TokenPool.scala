package cromwell.engine.workflow.tokens

import java.util.UUID

import cromwell.core.JobExecutionToken
import JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.TokenPool.TokenPoolPop

import scala.language.postfixOps

sealed trait TokenPool {
  def currentLoans: Set[JobExecutionToken]
  def pop(): TokenPoolPop
  def push(jobExecutionToken: JobExecutionToken): TokenPool
}

object TokenPool {

  case class TokenPoolPop(newTokenPool: TokenPool, poppedItem: Option[JobExecutionToken])

  def apply(tokenType: JobExecutionTokenType): TokenPool = {
    tokenType.maxPoolSize map { ps =>
      val pool = (1 to ps toList) map { _ => JobExecutionToken(tokenType, UUID.randomUUID()) }
      SizedTokenPool(pool.toSet, Set.empty)
    } getOrElse {
      InfiniteTokenPool(tokenType, Set.empty)
    }
  }
}

case class SizedTokenPool(pool: Set[JobExecutionToken], override val currentLoans: Set[JobExecutionToken]) extends TokenPool {

  override def pop(): TokenPoolPop = pool.toList match {
    case head :: tail => TokenPoolPop(SizedTokenPool(tail.toSet, currentLoans + head), Some(head))
    case Nil => TokenPoolPop(SizedTokenPool(Set.empty, currentLoans), None)
  }

  override def push(token: JobExecutionToken) = {
    if (currentLoans.contains(token)) { SizedTokenPool(pool + token, currentLoans - token) } else this
  }
}

case class InfiniteTokenPool(tokenType: JobExecutionTokenType, override val currentLoans: Set[JobExecutionToken]) extends TokenPool {
  override def pop() = {
    val newToken = JobExecutionToken(tokenType, UUID.randomUUID())
    TokenPoolPop(InfiniteTokenPool(tokenType, currentLoans + newToken), Some(newToken))
  }
  override def push(token: JobExecutionToken) = if (currentLoans.contains(token)) { InfiniteTokenPool(tokenType, currentLoans - token) } else this
}
