package cromwell.database.sql.tables

import java.sql.Timestamp

case class WorkflowStoreEntry
(
  workflowUuid: String,
  workflowSource: String,
  workflowInputs: String,
  workflowOptions: String,
  state: String,
  timestamp: Timestamp,
  workflowStoreId: Option[Int] = None
)
