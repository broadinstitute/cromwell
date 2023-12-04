package cromwell.database.sql.tables

case class SubWorkflowStoreEntry(
  rootWorkflowId: Option[Long],
  parentWorkflowExecutionUuid: String,
  callFullyQualifiedName: String,
  callIndex: Int,
  callAttempt: Int,
  subWorkflowExecutionUuid: String,
  subWorkflowStoreEntryId: Option[Long] = None
)
