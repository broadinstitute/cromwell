package wdl4s.wom.graph

import wdl4s.wdl.types.WdlType
import wdl4s.wom.graph.GraphNodePort.{ConnectedInputPort, OutputPort}

final case class GraphOutputNode(name: String, womType: WdlType, source: OutputPort) extends GraphNode {
  override def inputPorts: Set[GraphNodePort.InputPort] = Set(ConnectedInputPort(name, womType, source, _ => this))
  override def outputPorts: Set[GraphNodePort.OutputPort] = Set.empty
}
