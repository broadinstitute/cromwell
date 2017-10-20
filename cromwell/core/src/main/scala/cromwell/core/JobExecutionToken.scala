package cromwell.core

import java.util.UUID

import cromwell.core.JobExecutionToken.JobExecutionTokenType

case class JobExecutionToken(jobExecutionTokenType: JobExecutionTokenType, id: UUID)

object JobExecutionToken {
  case class JobExecutionTokenType(backend: String, maxPoolSize: Option[Int])
}
