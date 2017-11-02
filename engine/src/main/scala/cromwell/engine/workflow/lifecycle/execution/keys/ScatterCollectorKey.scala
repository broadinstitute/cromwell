package cromwell.engine.workflow.lifecycle.execution.keys

import common.validation.ErrorOr.ErrorOr
import cromwell.core.{ExecutionStatus, JobKey}
import cromwell.engine.workflow.lifecycle.execution.{WorkflowExecutionActorData, WorkflowExecutionDiff}
import wom.graph.GraphNodePort.ScatterGathererPort
import wom.graph.{GraphNode, ScatterNode}

/**
  * Key that becomes runnable when all shards of a collectible node are complete and need to be collected to form the output of this
  * call outside the scatter block.
  */
private [execution] case class ScatterCollectorKey(node: GraphNode,
                                                   scatterGatherPorts: Set[ScatterGathererPort],
                                                   scatter: ScatterNode,
                                                   scatterWidth: Int) extends JobKey {
  override val index = None
  override val attempt = 1
  override val tag = s"Collector-${node.localName}"

  def processRunnable(data: WorkflowExecutionActorData): ErrorOr[WorkflowExecutionDiff] = {
    data.valueStore.collectShards(this) map { outputs =>
      // TODO WOM: https://github.com/broadinstitute/cromwell/issues/2616
      //        val adjustedOutputs: CallOutputs = if (isInBypassed) {
      //          outputs map {
      //            case (outputKey, value) => outputKey -> WomOptionalValue.none(output._2.womValue.womType
      //          }
      //        } else outputs
      WorkflowExecutionDiff(
        executionStoreChanges = Map(this -> ExecutionStatus.Done),
        valueStoreAdditions = outputs
      )
    }
  }
}
