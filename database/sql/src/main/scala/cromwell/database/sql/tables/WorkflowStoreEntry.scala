package cromwell.database.sql.tables

import java.sql.{Blob, Clob, Timestamp}

import cromwell.database.sql.tables.WorkflowStoreEntry.WorkflowStoreState.WorkflowStoreState

object WorkflowStoreEntry {
  object WorkflowStoreState extends Enumeration {
    type WorkflowStoreState = Value
    val Submitted = Value("Submitted")
    val Running = Value("Running")
    val Aborting = Value("Aborting")
    val OnHold = Value("On Hold")
  }
}

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
  workflowState: WorkflowStoreState,
  submissionTime: Timestamp,
  importsZip: Option[Blob],
  customLabels: Clob,
  cromwellId: Option[String],
  heartbeatTimestamp: Option[Timestamp],
  workflowStoreEntryId: Option[Int] = None
)
