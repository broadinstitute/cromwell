package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import cromwell.core.{HogGroup, WorkflowId, WorkflowSourceFilesCollection}

/**
  * States of a workflow for which it can be fetched from the workflow store and started.
  */
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

final case class WorkflowToStart(id: WorkflowId,
                                 submissionTime: OffsetDateTime,
                                 sources: WorkflowSourceFilesCollection,
                                 state: StartableState,
                                 hogGroup: HogGroup)
