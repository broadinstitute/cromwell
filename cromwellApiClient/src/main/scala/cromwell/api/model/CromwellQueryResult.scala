package cromwell.api.model

import java.time.OffsetDateTime
import spray.json.DefaultJsonProtocol
import cromwell.api.model.WorkflowIdJsonFormatter._
import cromwell.api.model.WorkflowStatusJsonFormatter._

case class CromwellQueryResults(results: Seq[CromwellQueryResult])

case class CromwellQueryResult(name: Option[String], id: WorkflowId, status: WorkflowStatus, end: Option[OffsetDateTime], start: Option[OffsetDateTime], metadataArchiveStatus: String)

object CromwellQueryResultJsonSupport extends DefaultJsonProtocol {
  implicit val CromwellQueryResultJsonFormat = jsonFormat6(CromwellQueryResult)
  implicit val CromwellQueryResultsJsonFormat = jsonFormat1(CromwellQueryResults)
}
