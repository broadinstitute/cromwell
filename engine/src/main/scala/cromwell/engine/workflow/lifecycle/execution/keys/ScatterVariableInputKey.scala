package cromwell.engine.workflow.lifecycle.execution.keys

import cromwell.core.JobKey
import wom.graph.GraphInputNode
import wom.values.WdlArray.WdlArrayLike

private [execution] case class ScatterVariableInputKey(node: GraphInputNode, wdlArrayLike: WdlArrayLike) extends JobKey {
  override def index: Option[Int] = None
  override def attempt: Int = 1
  override def tag: String = node.localName
}
