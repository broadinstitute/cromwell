package wom.graph

import cats.instances.list._
import cats.syntax.apply._
import cats.syntax.functor._
import cats.syntax.traverse._
import cats.syntax.validated._
import lenthall.collections.EnhancedCollections._
import lenthall.validation.ErrorOr.ErrorOr
import shapeless.{:+:, CNil}
import wdl.values.WdlValue
import wom.expression.WomExpression
import wom.graph.GraphNodePort.InputPort

/**
  * A sealed set of graph nodes.
  */
final case class Graph private(nodes: Set[GraphNode]) {
  lazy val inputNodes: Set[GraphInputNode] = nodes.filterByType[GraphInputNode]
  lazy val externalInputNodes: Set[ExternalGraphInputNode] = nodes.filterByType[ExternalGraphInputNode]
  lazy val outputNodes: Set[GraphOutputNode] = nodes.filterByType[GraphOutputNode]
  lazy val calls: Set[CallNode] = nodes.filterByType[CallNode]
  lazy val scatters: Set[ScatterNode] = nodes.filterByType[ScatterNode]

  def outputByName(name: String): Option[GraphOutputNode] = outputNodes.find(_.name == name)
}

object Graph {

  type ResolvedExecutableInput = WdlValue :+: WomExpression :+: CNil

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
