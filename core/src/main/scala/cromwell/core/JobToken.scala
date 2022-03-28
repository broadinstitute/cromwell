package cromwell.core

import java.util.UUID
import _root_.io.circe.generic.JsonCodec
import cromwell.core.JobToken.JobTokenType

case class JobToken(jobTokenType: JobTokenType, id: UUID)

object JobToken {
  @JsonCodec
  case class JobTokenType(backend: String, maxPoolSize: Option[Int], hogFactor: Int)
}
