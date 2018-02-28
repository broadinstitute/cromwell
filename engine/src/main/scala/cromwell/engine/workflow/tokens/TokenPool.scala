package cromwell.engine.workflow.tokens

import java.util.UUID

import cromwell.core.JobExecutionToken
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import io.github.andrebeat.pool._

object TokenPool {
  // 10 Million
  val MaxCapacity = 10000000
}

final case class TokenPool(tokenType: JobExecutionTokenType) extends SimplePool[JobExecutionToken](
  capacity = tokenType.maxPoolSize.getOrElse(TokenPool.MaxCapacity),
  referenceType = ReferenceType.Strong,
  _factory = () => JobExecutionToken(tokenType, UUID.randomUUID()),
  _reset = Function.const(()),
  _dispose = Function.const(()),
  _healthCheck = Function.const(true)
)
