package cromwell.engine.workflow.lifecycle.execution.keys

import cromwell.core.JobKey
import wom.graph.GraphNode

private [execution] case class WorkflowInputKey(node: GraphNode) extends JobKey {
  override def index: Option[Int] = None
  override def attempt: Int = 1
  override def tag: String = node.localName
}
