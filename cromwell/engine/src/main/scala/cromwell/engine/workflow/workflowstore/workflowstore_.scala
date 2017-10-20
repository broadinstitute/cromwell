package cromwell.engine.workflow.workflowstore

import cromwell.core.{WorkflowId, WorkflowSourceFilesCollection}
import cromwell.engine.workflow.workflowstore.WorkflowStoreState.StartableState

sealed trait WorkflowStoreState {def isStartable: Boolean}

object WorkflowStoreState {
  case object Running extends WorkflowStoreState { override def isStartable = false }
  sealed trait StartableState extends WorkflowStoreState { override def isStartable = true }
  case object Submitted extends StartableState
  case object Restartable extends StartableState
}

final case class WorkflowToStart(id: WorkflowId, sources: WorkflowSourceFilesCollection, state: StartableState)
