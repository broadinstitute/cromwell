package cromwell.engine.workflow

import cromwell.core.{WorkflowId, WorkflowSourceFiles}
import cromwell.engine.workflow.workflowstore.WorkflowStoreState.StartableState

package object workflowstore {

  sealed trait WorkflowStoreState {def isStartable: Boolean}

  object WorkflowStoreState {
    case object Running extends WorkflowStoreState { override def isStartable = false }
    sealed trait StartableState extends WorkflowStoreState { override def isStartable = true }
    case object Submitted extends StartableState
    case object Restartable extends StartableState
  }

  final case class WorkflowToStart(id: WorkflowId, sources: WorkflowSourceFiles, state: StartableState)
}
