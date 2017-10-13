package wom.graph

import wom.types.{WdlArrayType, WdlOptionalType, WdlType}

sealed trait GraphNodePort {
  def name: String
  def womType: WdlType

  /**
    * The GraphNode which owns this port.
    */
  def graphNode: GraphNode
}

object GraphNodePort {

  // TODO: It'd be really cool if these could be typed (eg InputPort[WdlString], OutputPort[WdlInteger] but
  // TODO: we'd have to think about coercion... maybe some sort of implicit CoercionSocket[WdlString, WdlInteger]...?
  sealed trait InputPort extends GraphNodePort {
    def upstream: OutputPort
  }

  sealed trait OutputPort extends GraphNodePort {
    def identifier: WomIdentifier

    /**
      * Alias for identifier.localName.asString
      */
    override def name = identifier.localName.value
    // TODO: Might end up wanting a backwards link to the InputPorts that use this (eg def downstream: Set[InputPort])?
  }

  /**
    * A mixin trait that allows a port with a 'g: Unit => GraphNode' be linked to its [[wom.graph.GraphNode]]
    * after the GraphNode is constructed, using the [[wom.graph.GraphNode.GraphNodeSetter]]
    */
  sealed trait DelayedGraphNodePort { this: GraphNodePort =>
    val g: Unit => GraphNode
    override lazy val graphNode = g.apply(())
  }

  final case class ConnectedInputPort(name: String, womType: WdlType, upstream: OutputPort, g: Unit => GraphNode) extends InputPort with DelayedGraphNodePort

  /**
    * For any graph node that uses a declarations to produce outputs (e.g. call, declaration):
    */
  object GraphNodeOutputPort {
    def apply(name: String, womType: WdlType, graphNode: GraphNode): GraphNodeOutputPort = {
      GraphNodeOutputPort(WomIdentifier(LocalName(name), graphNode.identifier.fullyQualifiedName.combine(name)), womType, graphNode)
    }
  }
  case class GraphNodeOutputPort(override val identifier: WomIdentifier, womType: WdlType, graphNode: GraphNode) extends OutputPort

  /**
    * Represents the gathered output from a call/declaration in a ScatterNode.
    */
  final case class ScatterGathererPort(womType: WdlArrayType, outputToGather: PortBasedGraphOutputNode, g: Unit => GraphNode) extends OutputPort with DelayedGraphNodePort {
    // Since this port just wraps a PortBasedGraphOutputNode which itself wraps an output port, we can re-use the same identifier
    override def identifier: WomIdentifier = outputToGather.identifier
  }

  /**
    * Represents the conditional output from a call or declaration in a ConditionalNode
    */
  final case class ConditionalOutputPort(womType: WdlOptionalType, outputToExpose: PortBasedGraphOutputNode, g: Unit => GraphNode) extends OutputPort with DelayedGraphNodePort {
    // Since this port just wraps a PortBasedGraphOutputNode which itself wraps an output port, we can re-use the same identifier
    override def identifier: WomIdentifier = outputToExpose.identifier
  }
}
