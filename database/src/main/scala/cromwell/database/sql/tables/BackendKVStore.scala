package cromwell.database.sql.tables

case class BackendKVStore
(
  workflowExecutionUuid: String,
  callFqn: String,
  callIndex: Int,
  callAttempt: Int,
  backendJobKey: String,
  backendJobValue: String,
  backendKVStoreID: Option[Int] = None
)
