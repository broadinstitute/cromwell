package cromwell.database.sql.tables

import java.sql.{Blob, Clob, Timestamp}

case class WorkflowStoreEntry
(
  workflowExecutionUuid: String,
  workflowDefinition: Clob,
  workflowInputs: Clob,
  workflowOptions: Clob,
  workflowState: String,
  submissionTime: Timestamp,
  importsZip: Option[Blob],
  customLabels: Clob,
  workflowStoreEntryId: Option[Int] = None
)
