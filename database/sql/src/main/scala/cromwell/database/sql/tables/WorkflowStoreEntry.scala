package cromwell.database.sql.tables

import java.sql.Timestamp

import javax.sql.rowset.serial.{SerialBlob, SerialClob}

case class WorkflowStoreEntry
(
  workflowExecutionUuid: String,
  workflowDefinition: Option[SerialClob],
  workflowUrl: Option[String],
  workflowRoot: Option[String],
  workflowType: Option[String],
  workflowTypeVersion: Option[String],
  workflowInputs: Option[SerialClob],
  workflowOptions: Option[SerialClob],
  workflowState: String,
  submissionTime: Timestamp,
  importsZip: Option[SerialBlob],
  customLabels: SerialClob,
  cromwellId: Option[String],
  heartbeatTimestamp: Option[Timestamp],
  hogGroup: Option[String],
  workflowStoreEntryId: Option[Int] = None
)
