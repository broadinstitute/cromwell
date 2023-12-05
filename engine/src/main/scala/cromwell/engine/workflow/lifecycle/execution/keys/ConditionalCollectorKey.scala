package cromwell.engine.workflow.lifecycle.execution.keys

import common.validation.ErrorOr.ErrorOr
import cromwell.core.ExecutionIndex.ExecutionIndex
import cromwell.core.{ExecutionStatus, JobKey}
import cromwell.engine.workflow.lifecycle.execution.{WorkflowExecutionActorData, WorkflowExecutionDiff}
import wom.graph.{GraphNode, PortBasedGraphOutputNode}
import wom.graph.GraphNodePort.ConditionalOutputPort

/**
  * Key that becomes runnable when a node inside a conditional node is complete.
  * This is needed so that the ConditionalOutputPort of the conditional can be given a value.
  */
private[execution] case class ConditionalCollectorKey(conditionalOutputPort: ConditionalOutputPort,
                                                      index: ExecutionIndex
) extends JobKey {
  val outputNodeToCollect: PortBasedGraphOutputNode = conditionalOutputPort.outputToExpose
  override val node: GraphNode = conditionalOutputPort.outputToExpose
  override val attempt = 1
  override val tag = s"Collector-${node.localName}"

  def processRunnable(data: WorkflowExecutionActorData): ErrorOr[WorkflowExecutionDiff] =
    data.valueStore.collectConditional(this) map { outputs =>
      WorkflowExecutionDiff(
        executionStoreChanges = Map(this -> ExecutionStatus.Done),
        valueStoreAdditions = outputs
      )
    }
}
