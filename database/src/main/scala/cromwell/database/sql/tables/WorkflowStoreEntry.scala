package cromwell.database.sql.tables

import java.sql.Timestamp

case class WorkflowStoreEntry
(
  workflowUuid: String,
  workflowSource: String,
  workflowInputs: Option[String],
  workflowOptions: Option[String],
  state: String, // TODONE: This is a database column, and **the column** has no restrictions on the string stored.
  timestamp: Timestamp,
  workflowStoreId: Option[Int] = None
)
