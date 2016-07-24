package cromwell.database.sql.tables

import java.sql.Timestamp

case class WorkflowMetadataSummary
(
  workflowUuid: String,
  name: Option[String],
  status: Option[String],
  startDate: Option[Timestamp],
  endDate: Option[Timestamp],
  workflowMetadataSummaryId: Option[Long] = None
)
