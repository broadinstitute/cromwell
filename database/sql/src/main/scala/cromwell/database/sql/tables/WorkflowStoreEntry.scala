package cromwell.database.sql.tables

import java.sql.{Clob, Timestamp}

case class WorkflowStoreEntry
(
  workflowExecutionUuid: String,
  workflowDefinition: Clob,
  workflowInputs: Clob,
  workflowOptions: Clob,
  workflowState: String,
  submissionTime: Timestamp,
  workflowStoreEntryId: Option[Int] = None
)
