package cromwell.database.sql.tables

import java.sql.Timestamp

case class WorkflowMetadataSummaryEntry
(
  workflowExecutionUuid: String,
  workflowName: Option[String],
  workflowStatus: Option[String],
  startTimestamp: Option[Timestamp],
  endTimestamp: Option[Timestamp],
  workflowMetadataSummaryEntryId: Option[Long] = None
)
