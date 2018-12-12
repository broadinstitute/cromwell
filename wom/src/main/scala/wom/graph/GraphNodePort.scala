package wom.graph

import wom.expression.WomExpression
import wom.types.{WomArrayType, WomNothingType, WomOptionalType, WomType}

sealed trait GraphNodePort {
  def name: String
  def womType: WomType

  /**
    * The GraphNode which owns this port.
    */
  def graphNode: GraphNode
}

object GraphNodePort {

  // TODO: It'd be really cool if these could be typed (eg InputPort[WomString], OutputPort[WomInteger] but
  // TODO: we'd have to think about coercion... maybe some sort of implicit CoercionSocket[WomString, WomInteger]...?
  sealed trait InputPort extends GraphNodePort {
    def upstream: OutputPort
  }

  sealed trait OutputPort extends GraphNodePort {
    def identifier: WomIdentifier
    override def toString: String = s"${getClass.getSimpleName}(${identifier.fullyQualifiedName.value})"

    /**
      * Alias for identifier.localName.asString
      */
    override def name = identifier.localName.value

    /**
      * The name with which to refer to this port from within the same Node (as opposed to 'name' which is how it is
      * referred to from the outside. Currently this is only different in WDL, so we can hard code this function:
      */
    def internalName = identifier.localName.value.split("\\.").last
  }

  /**
    * A mixin trait that allows a port with a 'g: Unit => GraphNode' be linked to its [[wom.graph.GraphNode]]
    * after the GraphNode is constructed, using the [[wom.graph.GraphNode.GraphNodeSetter]]
    */
  sealed trait DelayedGraphNodePort { this: GraphNodePort =>
    val g: Unit => GraphNode
    override lazy val graphNode = g.apply(())
  }

  final case class ConnectedInputPort(name: String, womType: WomType, upstream: OutputPort, g: Unit => GraphNode) extends InputPort with DelayedGraphNodePort

  /**
    * For any graph node that uses a declarations to produce outputs (e.g. call, declaration):
    */
  object GraphNodeOutputPort {
    def apply(name: String, womType: WomType, graphNode: GraphNode): GraphNodeOutputPort = {
      GraphNodeOutputPort(WomIdentifier(LocalName(name), graphNode.identifier.fullyQualifiedName.combine(name)), womType, graphNode)
    }
  }
  case class GraphNodeOutputPort(override val identifier: WomIdentifier, womType: WomType, graphNode: GraphNode) extends OutputPort

  object ExpressionBasedOutputPort
  case class ExpressionBasedOutputPort(override val identifier: WomIdentifier, womType: WomType, graphNode: GraphNode, expression: WomExpression) extends OutputPort

  /**
    * Represents the gathered output from a call/declaration in a ScatterNode.
    */
  final case class ScatterGathererPort(womType: WomArrayType, outputToGather: PortBasedGraphOutputNode, g: Unit => GraphNode) extends OutputPort with DelayedGraphNodePort {
    // Since this port just wraps a PortBasedGraphOutputNode which itself wraps an output port, we can re-use the same identifier
    override def identifier: WomIdentifier = outputToGather.identifier
  }

  /**
    * Represents the conditional output from a call or declaration in a ConditionalNode
    */
  final case class ConditionalOutputPort(outputToExpose: PortBasedGraphOutputNode, g: Unit => ConditionalNode) extends OutputPort with DelayedGraphNodePort {
    // Since this port just wraps a PortBasedGraphOutputNode which itself wraps an output port, we can re-use the same identifier
    override def identifier: WomIdentifier = outputToExpose.identifier
    override val womType: WomOptionalType = WomOptionalType(outputToExpose.womType).flatOptionalType
    lazy val conditionalNode: ConditionalNode = g(())
  }

  final case class NodeCompletionPort(g: Unit => GraphNode) extends OutputPort with DelayedGraphNodePort {
    override lazy val identifier: WomIdentifier = {
      val name = "__after"
      WomIdentifier(LocalName(name), graphNode.identifier.fullyQualifiedName.combine(name))
    }
    override def womType: WomType = WomNothingType
  }

  /**
    * Represents an output port from a workflow call, based on the output that it exposes.
    */
  final case class SubworkflowCallOutputPort(identifier: WomIdentifier, outputToExpose: GraphOutputNode, workflowCallNode: WorkflowCallNode) extends OutputPort {
    override val womType: WomType = outputToExpose.womType
    override val graphNode: GraphNode = workflowCallNode
  }
}
