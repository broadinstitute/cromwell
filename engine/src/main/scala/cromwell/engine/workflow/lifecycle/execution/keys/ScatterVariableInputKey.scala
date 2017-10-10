package cromwell.engine.workflow.lifecycle.execution.keys

import cromwell.core.JobKey
import wdl.values.WdlArray.WdlArrayLike
import wom.graph.GraphInputNode

private [execution] case class ScatterVariableInputKey(node: GraphInputNode, wdlArrayLike: WdlArrayLike) extends JobKey {
  override def index: Option[Int] = None
  override def attempt: Int = 1
  override def tag: String = node.localName
}
