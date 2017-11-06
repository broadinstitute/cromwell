package cromwell.database.sql.tables

import java.sql.{Blob, Clob, Timestamp}

import cromwell.database.sql.tables.WorkflowStoreEntry.WorkflowStoreState.WorkflowStoreState

object WorkflowStoreEntry {
  type StateAndRestarted = (String, Boolean)

  object WorkflowStoreState extends Enumeration {
    type WorkflowStoreState = Value
    val Submitted = Value("Submitted")
    val Running = Value("Running")
    val Aborting = Value("Aborting")
  }
}

case class WorkflowStoreEntry
(
  workflowExecutionUuid: String,
  workflowDefinition: Option[Clob],
  workflowType: Option[String],
  workflowTypeVersion: Option[String],
  workflowInputs: Option[Clob],
  workflowOptions: Option[Clob],
  workflowState: WorkflowStoreState,
  restarted: Boolean,
  submissionTime: Timestamp,
  importsZip: Option[Blob],
  customLabels: Clob,
  workflowStoreEntryId: Option[Int] = None
)
