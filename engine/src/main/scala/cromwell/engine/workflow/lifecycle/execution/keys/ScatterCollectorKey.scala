package cromwell.engine.workflow.lifecycle.execution.keys

import common.validation.ErrorOr.ErrorOr
import cromwell.core.{ExecutionStatus, JobKey}
import cromwell.engine.workflow.lifecycle.execution.{WorkflowExecutionActorData, WorkflowExecutionDiff}
import wom.graph.GraphNodePort.ScatterGathererPort
import wom.graph.ScatterNode.ScatterCollectionFunction
import wom.graph.{GraphNode, PortBasedGraphOutputNode}

/**
  * Key that becomes runnable when all shards of a collectible node are complete and need to be collected to form the output of this
  * call outside the scatter block.
  */
private[execution] case class ScatterCollectorKey(scatterGatherPort: ScatterGathererPort,
                                                  scatterWidth: Int,
                                                  scatterCollectionFunction: ScatterCollectionFunction
) extends JobKey {
  val outputNodeToGather: PortBasedGraphOutputNode = scatterGatherPort.outputToGather
  override val node: GraphNode = outputNodeToGather
  override val index = None
  override val attempt = 1
  override val tag = s"Collector-${node.localName}"

  def processRunnable(data: WorkflowExecutionActorData): ErrorOr[WorkflowExecutionDiff] =
    data.valueStore.collectShards(this) map { outputs =>
      WorkflowExecutionDiff(
        executionStoreChanges = Map(this -> ExecutionStatus.Done),
        valueStoreAdditions = outputs
      )
    }
}
