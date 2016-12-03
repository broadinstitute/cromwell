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
  importsZipFile: Option[Blob],
  workflowStoreEntryId: Option[Int] = None
)
