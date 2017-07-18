package wdl4s.wom.graph

import wdl4s.wdl.types.{WdlOptionalType, WdlType}
import wdl4s.wom.expression.Expression
import wdl4s.wom.graph.GraphNodePort.GraphNodeOutputPort

sealed trait GraphInputNode extends GraphNode {
  def name: String
  def womType: WdlType
  val singleOutputPort: GraphNodeOutputPort = GraphNodeOutputPort(name, womType, this)

  override val inputPorts: Set[GraphNodePort.InputPort] = Set.empty
  override val outputPorts: Set[GraphNodePort.OutputPort] = Set(singleOutputPort)
}

final case class RequiredGraphInputNode(name: String, womType: WdlType) extends GraphInputNode
final case class OptionalGraphInputNode(name: String, womType: WdlOptionalType) extends GraphInputNode
final case class OptionalGraphInputNodeWithDefault(name: String, womType: WdlType, default: Expression) extends GraphInputNode
