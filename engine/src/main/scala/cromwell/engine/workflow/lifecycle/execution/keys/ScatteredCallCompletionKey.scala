package cromwell.engine.workflow.lifecycle.execution.keys

import cats.syntax.traverse._
import cats.instances.list._
import common.validation.ErrorOr.ErrorOr
import cromwell.core.{ExecutionStatus, JobKey}
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore.ValueKey
import cromwell.engine.workflow.lifecycle.execution.{WorkflowExecutionActorData, WorkflowExecutionDiff}
import wom.graph.{CallNode, GraphNode}
import wom.values.{WomArray, WomValue}

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
    val outputLookup: ErrorOr[List[WomValue]] = (0 until scatterWidth).toList.traverse(index => data.valueStore.resolve(Option(index))(call.singleCompletionPort))

    outputLookup map { outputList =>
      val r = WomArray(outputList)
      WorkflowExecutionDiff(
        executionStoreChanges = Map(this -> ExecutionStatus.Done),
        valueStoreAdditions = Map(ValueKey(call.singleCompletionPort, None) -> r)
      )
    }
  }
}
