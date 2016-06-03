package cromwell.database.obj

import java.sql.Timestamp

case class WorkflowMetadataSummary
(
  workflowUuid: String,
  name: Option[String] = None,
  status: Option[String] = None,
  startDate: Option[Timestamp] = None,
  endDate: Option[Timestamp] = None,
  workflowMetadataSummaryId: Option[Long] = None
)
