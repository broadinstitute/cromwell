package cromwell.database.sql.tables

case class CallCachingEntry
(
  workflowExecutionUuid: String,
  callFullyQualifiedName: String,
  jobIndex: Int,
  jobAttempt: Option[Int],
  returnCode: Option[Int],
  allowResultReuse: Boolean,
  callCachingEntryId: Option[Int] = None
)
