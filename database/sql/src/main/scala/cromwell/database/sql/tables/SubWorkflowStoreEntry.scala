package cromwell.database.sql.tables

case class SubWorkflowStoreEntry
(
  rootWorkflowId: Option[Int],
  parentWorkflowExecutionUuid: String,
  callFullyQualifiedName: String,
  callIndex: Int,
  callAttempt: Int,
  subWorkflowExecutionUuid: String,
  subWorkflowStoreEntryId: Option[Int] = None
)
