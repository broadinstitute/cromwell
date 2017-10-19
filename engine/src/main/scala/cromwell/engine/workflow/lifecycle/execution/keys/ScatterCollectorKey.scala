package cromwell.engine.workflow.lifecycle.execution.keys

import cromwell.core.JobKey
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
}
