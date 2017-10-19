package cromwell.engine.workflow.lifecycle.execution.keys

import cromwell.core.ExecutionIndex.ExecutionIndex
import cromwell.core.JobKey
import wom.graph.ConditionalNode
import wom.graph.GraphNodePort.ConditionalOutputPort

/**
  * Key that becomes runnable when a node inside a conditional node is complete.
  * This is needed so that the ConditionalOutputPort of the conditional can be given a value.
  */
private [execution] case class ConditionalCollectorKey(conditionalOutputPort: ConditionalOutputPort,
                                                       index: ExecutionIndex,
                                                       conditional: ConditionalNode) extends JobKey {
  override val node = conditionalOutputPort.outputToExpose
  override val attempt = 1
  override val tag = s"Collector-${node.localName}"
}
