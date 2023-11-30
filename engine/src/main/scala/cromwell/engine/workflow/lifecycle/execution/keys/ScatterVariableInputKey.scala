package cromwell.engine.workflow.lifecycle.execution.keys

import cromwell.core.JobKey
import wom.graph.GraphInputNode
import wom.values.WomArray.WomArrayLike

private[execution] case class ScatterVariableInputKey(node: GraphInputNode, womArrayLike: WomArrayLike) extends JobKey {
  override def index: Option[Int] = None
  override def attempt: Int = 1
  override def tag: String = node.localName
}
