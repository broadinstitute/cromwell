package cromwell.database.sql.tables

import java.sql.{Blob, Clob, Timestamp}

case class WorkflowStoreEntry
(
  workflowExecutionUuid: String,
  workflowDefinition: Option[Clob],
  workflowUrl: Option[String],
  workflowRoot: Option[String],
  workflowType: Option[String],
  workflowTypeVersion: Option[String],
  workflowInputs: Option[Clob],
  workflowOptions: Option[Clob],
  workflowState: String,
  submissionTime: Timestamp,
  importsZip: Option[Blob],
  customLabels: Clob,
  cromwellId: Option[String],
  heartbeatTimestamp: Option[Timestamp],
  workflowStoreEntryId: Option[Int] = None
)
