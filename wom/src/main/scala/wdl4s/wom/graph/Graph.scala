package wdl4s.wom.graph

import cats.syntax.validated._
import cats.syntax.cartesian._
import cats.instances.list._
import cats.syntax.functor._
import cats.syntax.traverse._

import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.wom.graph.GraphNodePort.InputPort

/**
  * A sealed set of graph nodes.
  */
final case class Graph private(nodes: Set[GraphNode]) {
  // FIXME: this actually seems pretty WDL specific, and even then specific to the case of a WDL that doesn't have
  // an explicit outputs block.  I'd remove this but the WDL/WOM test relies on it.
  def withDefaultOutputs: Graph = {
    val defaultOutputs: Set[GraphNode] = (nodes collect {
      case callNode: CallNode => callNode.outputPorts.map(op => PortBasedGraphOutputNode(s"${callNode.name}.${op.name}", op.womType, op))
    }).flatten
    Graph(nodes.union(defaultOutputs))
  }
}

object Graph {

  /**
    * Checks that every input port for every node in the graph references an upstream node that is also in the graph.
    * Assuming it validates, construct the Graph case class.
    */
  def validateAndConstruct(nodes: Set[GraphNode]): ErrorOr[Graph] = {

    def boolToErrorOr(bool: Boolean, msg: => String): ErrorOr[Unit] = if (bool) ().validNel else msg.invalidNel

    def upstreamNodeInGraph(port: InputPort): ErrorOr[Unit] = {
      val upstreamOutputPort = port.upstream
      boolToErrorOr(nodes.exists(_ eq upstreamOutputPort.graphNode), s"The input link ${port.name} is linked to a node outside the graph set (${upstreamOutputPort.name})")
    }

    def portProperlyEmbedded(port: GraphNodePort, portFinder: GraphNode => Set[_ <: GraphNodePort]): ErrorOr[Unit] = {
      boolToErrorOr(portFinder(port.graphNode).exists(_ eq port), s"The port ${port.name} thinks it belongs to a Node (${port.graphNode}), but that Node doesn't think it owns it.")
    }

    def goodLink(port: InputPort): ErrorOr[Unit] = {
      val upstreamNodeValidation = upstreamNodeInGraph(port)
      val inputPortEmbeddedValidation = portProperlyEmbedded(port, _.inputPorts)
      val upstreamPortEmbeddedValidation = portProperlyEmbedded(port.upstream, _.outputPorts)

      (upstreamNodeValidation |@| inputPortEmbeddedValidation |@| upstreamPortEmbeddedValidation).tupled.void
    }

    def validateNode(node: GraphNode): ErrorOr[Unit] = {
      node.inputPorts.toList.traverse(goodLink).void
    }

    nodes.toList.traverse(validateNode).map(_ => Graph(nodes))
  }
}
