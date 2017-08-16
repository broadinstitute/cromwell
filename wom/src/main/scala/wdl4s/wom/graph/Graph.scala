package wdl4s.wom.graph

import cats.implicits._
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
    def goodLink(port: InputPort): ErrorOr[Unit] = {
      val upstreamOutputPort = port.upstream

      if (nodes.contains(upstreamOutputPort.graphNode)) {
        ().validNel
      } else {
        s"The input link ${port.name} is linked to a node outside the graph set (${upstreamOutputPort.name})".invalidNel
      }
    }

    def validateNode(node: GraphNode): ErrorOr[Unit] = {
      node.inputPorts.toList.traverse(goodLink).void
    }

    nodes.toList.traverse(validateNode).map(_ => Graph(nodes))
  }
}
