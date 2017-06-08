package cromwell.database.sql.tables

import java.sql.{Blob, Clob, Timestamp}

case class WorkflowStoreEntry
(
  workflowExecutionUuid: String,
  workflowDefinition: Option[Clob],
  workflowType: String,
  workflowTypeVersion: Option[String],
  workflowInputs: Option[Clob],
  workflowOptions: Option[Clob],
  workflowState: String,
  submissionTime: Timestamp,
  importsZip: Option[Blob],
  customLabels: Clob,
  workflowStoreEntryId: Option[Int] = None
)
