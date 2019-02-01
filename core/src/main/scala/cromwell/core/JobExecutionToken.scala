package cromwell.core

import java.util.UUID
import _root_.io.circe.generic.JsonCodec
import cromwell.core.JobExecutionToken.JobExecutionTokenType

case class JobExecutionToken(jobExecutionTokenType: JobExecutionTokenType, id: UUID)

object JobExecutionToken {
  @JsonCodec
  case class JobExecutionTokenType(backend: String, maxPoolSize: Option[Int], hogFactor: Int)
}
