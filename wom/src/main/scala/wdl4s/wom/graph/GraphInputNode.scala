package wdl4s.wom.graph

import wdl4s.wdl.types.{WdlOptionalType, WdlType}
import wdl4s.wom.expression.WomExpression
import wdl4s.wom.graph.GraphNodePort.GraphNodeOutputPort

sealed trait GraphInputNode extends GraphNode {
  def womType: WdlType
  val singleOutputPort: GraphNodeOutputPort = GraphNodeOutputPort(name, womType, this)

  override val inputPorts: Set[GraphNodePort.InputPort] = Set.empty
  override val outputPorts: Set[GraphNodePort.OutputPort] = Set(singleOutputPort)
}

final case class RequiredGraphInputNode(override val name: String, womType: WdlType) extends GraphInputNode
final case class OptionalGraphInputNode(override val name: String, womType: WdlOptionalType) extends GraphInputNode
final case class OptionalGraphInputNodeWithDefault(override val name: String, womType: WdlType, default: WomExpression) extends GraphInputNode
