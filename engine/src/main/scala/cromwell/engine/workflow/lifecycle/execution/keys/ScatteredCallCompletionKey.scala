package cromwell.engine.workflow.lifecycle.execution.keys

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import cromwell.core.{ExecutionStatus, JobKey}
import cromwell.engine.workflow.lifecycle.execution.{WorkflowExecutionActorData, WorkflowExecutionDiff}
import wom.graph.{CallNode, GraphNode}

/**
  * Key that should become runnable when all shards of a scattered call are complete.
  */
private [execution] case class ScatteredCallCompletionKey(call: CallNode,
                                                          scatterWidth: Int) extends JobKey {
  override val node: GraphNode = call
  override val index = None
  override val attempt = 1
  override val totalIndices = scatterWidth
  override val tag = s"CallCompletion-${node.localName}"

  def processRunnable(data: WorkflowExecutionActorData): ErrorOr[WorkflowExecutionDiff] = {
    WorkflowExecutionDiff(executionStoreChanges = Map(this -> ExecutionStatus.Done)).validNel
  }
}
