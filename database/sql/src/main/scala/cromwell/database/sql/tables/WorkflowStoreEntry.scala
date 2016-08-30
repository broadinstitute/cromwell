package cromwell.database.sql.tables

import java.sql.{Clob, Timestamp}

case class WorkflowStoreEntry
(
  workflowUuid: String,
  workflowSource: Clob,
  workflowInputs: Clob,
  workflowOptions: Clob,
  state: String,
  timestamp: Timestamp,
  workflowStoreId: Option[Int] = None
)
