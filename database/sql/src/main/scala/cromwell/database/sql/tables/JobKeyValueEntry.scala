package cromwell.database.sql.tables

case class BackendKVStoreEntry
(
  workflowExecutionUuid: String,
  callFqn: String,
  jobScatterIndex: Int,
  jobRetryAttempt: Int,
  storeKey: String,
  storeValue: String,
  backendKVStoreID: Option[Int] = None
)
