package cromwell.api.model

import java.time.OffsetDateTime
import spray.json.DefaultJsonProtocol
import cromwell.api.model.WorkflowIdJsonFormatter._
import cromwell.api.model.WorkflowStatusJsonFormatter._

case class CromwellQueryResults(results: Seq[CromwellQueryResult])

case class CromwellQueryResult(name: String, id: WorkflowId, status: WorkflowStatus, end: OffsetDateTime, start: OffsetDateTime)

object CromwellQueryResultJsonFormatter extends DefaultJsonProtocol {
  implicit val CromwellQueryResultJsonFormat = jsonFormat5(CromwellQueryResult)
  implicit val CromwellQueryResultsJsonFormat = jsonFormat1(CromwellQueryResults)
}
