package cromwell.engine.workflow.lifecycle.execution

import cromwell.core.ExecutionIndex._
import wdl4s.wdl.types.WdlType
import wdl4s.wdl.values.WdlValue
import wdl4s.wom.graph.GraphNode
import wdl4s.wom.graph.GraphNodePort.OutputPort

object WomOutputStore {
  case class OutputEntry(name: String, wdlType: WdlType, wdlValue: Option[WdlValue])
  case class OutputCallKey(call: GraphNode, index: ExecutionIndex)
  def empty = OutputStore(Map.empty)
}

case class WomOutputStore(store: Map[OutputPort, WdlValue]) {

  override def toString = store.map { case (k, l) => s"$k -> $l" } mkString System.lineSeparator

  def add(values: Map[OutputPort, WdlValue]) = this.copy(store = store ++ values)

  def get(outputPort: OutputPort): Option[WdlValue] = store.get(outputPort)
}
