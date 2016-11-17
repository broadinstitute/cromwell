package cromwell.database.sql.tables

case class SubWorkflowStoreEntry
(
  rootWorkflowExecutionUuid: String,
  parentWorkflowExecutionUuid: String,
  callFullyQualifiedName: String,
  jobIndex: Int,
  jobAttempt: Int,
  subWorkflowExecutionUuid: String,
  subWorkflowStoreEntryId: Option[Int] = None
)
