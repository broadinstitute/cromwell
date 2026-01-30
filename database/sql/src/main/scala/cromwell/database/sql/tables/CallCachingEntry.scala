package cromwell.database.sql.tables

import java.sql.Timestamp

case class CallCachingEntry(
  workflowExecutionUuid: String,
  callFullyQualifiedName: String,
  jobIndex: Int,
  jobAttempt: Option[Int],
  returnCode: Option[Int],
  allowResultReuse: Boolean,
  callCachingEntryId: Option[Long] = None,
  createdAt: Timestamp
)
