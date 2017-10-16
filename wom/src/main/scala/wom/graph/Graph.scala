package wom.graph

import cats.data.NonEmptyList
import cats.instances.list._
import cats.syntax.apply._
import cats.syntax.functor._
import cats.syntax.traverse._
import cats.syntax.validated._
import lenthall.collections.EnhancedCollections._
import lenthall.validation.ErrorOr.ErrorOr
import shapeless.{:+:, CNil}
import wom.expression.WomExpression
import wom.graph.GraphNodePort.InputPort
import wom.values.WomValue

/**
  * A sealed set of graph nodes.
  */
final case class Graph private(nodes: Set[GraphNode]) {
  lazy val inputNodes: Set[GraphInputNode] = nodes.filterByType[GraphInputNode]
  lazy val externalInputNodes: Set[ExternalGraphInputNode] = nodes.filterByType[ExternalGraphInputNode]
  lazy val outputNodes: Set[GraphOutputNode] = nodes.filterByType[GraphOutputNode]
  lazy val calls: Set[CallNode] = nodes.filterByType[CallNode]
  lazy val scatters: Set[ScatterNode] = nodes.filterByType[ScatterNode]

  def outputByName(name: String): Option[GraphOutputNode] = outputNodes.find(_.localName == name)
}

object Graph {

  type ResolvedExecutableInput = WomValue :+: WomExpression :+: CNil

  /**
    * Checks that every input port for every node in the graph references an upstream node that is also in the graph.
    * Assuming it validates, construct the Graph case class.
    */
  def validateAndConstruct(nodes: Set[GraphNode]): ErrorOr[Graph] = {

    def boolToErrorOr(bool: Boolean, msg: => String): ErrorOr[Unit] = if (bool) ().validNel else msg.invalidNel

    def upstreamNodeInGraph(port: InputPort): ErrorOr[Unit] = {
      val upstreamOutputPort = port.upstream
      boolToErrorOr(nodes.exists(_ eq upstreamOutputPort.graphNode), s"The input link ${port.name} on ${port.graphNode.localName} is linked to a node outside the graph set (${upstreamOutputPort.name})")
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

    // from https://stackoverflow.com/a/24729587/1498572
    def fqnUniqueness: ErrorOr[Unit] = nodes
      .collect({
        case callNode: CallNode => callNode
        case gin: GraphInputNode => gin
        case gon: GraphOutputNode => gon
      })  
      .toList // Important since nodes is a Set, we don't want duplicates to disappear automatically when mapping to FQN
      .map(_.identifier.fullyQualifiedName)
      .groupBy(identity)
      .collect({
        case (fqn, list) if list.lengthCompare(1) > 0 => fqn
      }).toList match {
      case Nil => ().validNel
      case head :: tail =>
        NonEmptyList.of(head, tail: _*).map(fqn => s"Two or more nodes have the same FullyQualifiedName: ${fqn.value}").invalid
    }

    (fqnUniqueness, nodes.toList.traverse(validateNode)) mapN { case (_, _) =>
      Graph(nodes)
    }
  }
}
