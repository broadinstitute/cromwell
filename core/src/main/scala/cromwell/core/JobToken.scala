package cromwell.core

import java.util.UUID
import spray.json.DefaultJsonProtocol._
import spray.json.RootJsonFormat
import cromwell.core.JobToken.JobTokenType

case class JobToken(jobTokenType: JobTokenType, id: UUID)

object JobToken {
  case class JobTokenType(backend: String, maxPoolSize: Option[Int], hogFactor: Int)

  object JobTokenType {
    implicit val jobTokenTypeFormat: RootJsonFormat[JobTokenType] = jsonFormat3(JobTokenType.apply)
  }
}
