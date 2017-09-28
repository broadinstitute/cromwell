package wom.graph

import cats.instances.list._
import cats.syntax.apply._
import cats.syntax.functor._
import cats.syntax.traverse._
import cats.syntax.validated._
import lenthall.collections.EnhancedCollections._
import lenthall.validation.ErrorOr.{ErrorOr, _}
import lenthall.validation.Validation._
import shapeless.{:+:, CNil, Coproduct}
import wdl.WorkflowRawInputs
import wdl.values.WdlValue
import wom.expression.WomExpression
import wom.graph.Graph.ResolvedWorkflowInput
import wom.graph.GraphNodePort.{InputPort, OutputPort}

/**
  * A sealed set of graph nodes.
  */
final case class Graph private(nodes: Set[GraphNode]) {
  lazy val inputNodes: Set[GraphInputNode] = nodes.filterByType[GraphInputNode]
  lazy val outputNodes: Set[GraphOutputNode] = nodes.filterByType[GraphOutputNode]
  lazy val calls: Set[CallNode] = nodes.filterByType[CallNode]
  lazy val scatters: Set[ScatterNode] = nodes.filterByType[ScatterNode]

  def outputByName(name: String): Option[GraphOutputNode] = outputNodes.find(_.name == name)

  /**
    * Maps GraphInputNode to their final value / expression using workflow inputs. Validate all required inputs are satisfied.
    * @param inputsMapping workflow inputs in the form of Map[FQN, Any]
    * @return validated mappings from GraphInputPort to ResolvedWorkflowInput, which can be a WdlValue or an expression
    */
  def validateWorkflowInputs(inputsMapping: WorkflowRawInputs): ErrorOr[Map[OutputPort, ResolvedWorkflowInput]] = {

    def coerceRawValue(value: Any, gin: ExternalGraphInputNode): ErrorOr[WdlValue] = {
      gin.womType.coerceRawValue(value).toErrorOr
    }
    
    def fromInputMapping(gin: ExternalGraphInputNode): Option[ErrorOr[ResolvedWorkflowInput]] = {
      inputsMapping.get(s"${gin.name}").map(coerceRawValue(_, gin).map(Coproduct[ResolvedWorkflowInput](_)))
    }

    def fallBack(gin: ExternalGraphInputNode): ErrorOr[ResolvedWorkflowInput] = gin match {
      case required: RequiredGraphInputNode => s"Cannot find an input value for ${required.name}".invalidNel
      case optionalWithDefault: OptionalGraphInputNodeWithDefault => Coproduct[ResolvedWorkflowInput](optionalWithDefault.default).validNel
      case optional: OptionalGraphInputNode => Coproduct[ResolvedWorkflowInput](optional.womType.none: WdlValue).validNel
    }

    nodes.collect({
      case gin: ExternalGraphInputNode => 
        // The compiler needs the type ascription for some reason
        (gin.singleOutputPort: OutputPort) -> fromInputMapping(gin).getOrElse(fallBack(gin))
    }).toMap.sequence
  }
}

object Graph {

  type ResolvedWorkflowInput = WdlValue :+: WomExpression :+: CNil

  /**
    * Checks that every input port for every node in the graph references an upstream node that is also in the graph.
    * Assuming it validates, construct the Graph case class.
    */
  def validateAndConstruct(nodes: Set[GraphNode]): ErrorOr[Graph] = {

    def boolToErrorOr(bool: Boolean, msg: => String): ErrorOr[Unit] = if (bool) ().validNel else msg.invalidNel

    def upstreamNodeInGraph(port: InputPort): ErrorOr[Unit] = {
      val upstreamOutputPort = port.upstream
      boolToErrorOr(nodes.exists(_ eq upstreamOutputPort.graphNode), s"The input link ${port.name} on ${port.graphNode.name} is linked to a node outside the graph set (${upstreamOutputPort.name})")
    }

    def portProperlyEmbedded(port: GraphNodePort, portFinder: GraphNode => Set[_ <: GraphNodePort]): ErrorOr[Unit] = {
      boolToErrorOr(portFinder(port.graphNode).exists(_ eq port), s"The port ${port.name} thinks it belongs to a Node (${port.graphNode}), but that Node doesn't think it owns it.")
    }

    def goodLink(port: InputPort): ErrorOr[Unit] = {
      val upstreamNodeValidation = upstreamNodeInGraph(port)
      val inputPortEmbeddedValidation = portProperlyEmbedded(port, _.inputPorts)
      val upstreamPortEmbeddedValidation = portProperlyEmbedded(port.upstream, _.outputPorts)

      (upstreamNodeValidation, inputPortEmbeddedValidation, upstreamPortEmbeddedValidation).tupled.void
    }

    def validateNode(node: GraphNode): ErrorOr[Unit] = {
      node.inputPorts.toList.traverse(goodLink).void
    }

    nodes.toList.traverse(validateNode).map(_ => Graph(nodes))
  }
}
