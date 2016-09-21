package cromwell.database.sql.tables

case class JobKeyValueEntry
(
  workflowExecutionUuid: String,
  callFullyQualifiedName: String,
  jobIndex: Int,
  jobAttempt: Int,
  storeKey: String,
  storeValue: String,
  jobKeyValueEntryId: Option[Int] = None
)
