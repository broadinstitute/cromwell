package wom.graph

import shapeless.Coproduct
import wom.graph.CallNode.InputDefinitionPointer
import wom.graph.GraphNodePort.OutputPort

trait GraphNodeWithSingleOutputPort extends GraphNode {
  def singleOutputPort: OutputPort

  /**
    * Can be used to use this node as an InputDefinitionPointer
    */
  lazy val inputDefinitionPointer = Coproduct[InputDefinitionPointer](singleOutputPort: OutputPort)
}
