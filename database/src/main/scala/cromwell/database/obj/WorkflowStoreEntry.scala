package cromwell.database.obj

import java.sql.Timestamp

case class WorkflowStoreEntry
(
  workflowUuid: String,
  workflowSource: String,
  workflowInputs: Option[String],
  workflowOptions: Option[String],
  state: String, // TODO: This could probably be a WorkflowStoreEntryState and I could use TypedTypes. On the other hand... String...
  timestamp: Timestamp,
  workflowStoreId: Option[Int]
)

sealed trait WorkflowStoreEntryState

object WorkflowStoreEntryState {
  case object Running extends WorkflowStoreEntryState
  case object Restartable extends WorkflowStoreEntryState
  case object Submitted extends WorkflowStoreEntryState
}
