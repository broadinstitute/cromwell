package cromwell.engine.workflow.workflowstore

import cromwell.core.{WorkflowId, WorkflowSourceFilesCollection}

sealed trait StartableState {
  def restarted: Boolean
}

case object Submitted extends StartableState {
  override val restarted = false
}

case object RestartableRunning extends StartableState {
  override val restarted = true
}

case object RestartableAborting extends StartableState {
  override val restarted = true
}

final case class WorkflowToStart(id: WorkflowId, sources: WorkflowSourceFilesCollection, state: StartableState)
