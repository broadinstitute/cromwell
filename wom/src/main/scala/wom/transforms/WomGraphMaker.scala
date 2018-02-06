package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.graph.{Graph, GraphNode, GraphNodePort}

trait WomGraphMaker[A] {
  def toWomGraph(a: A, includeGraphNodes: Set[GraphNode], outerLookup: Map[String, GraphNodePort.OutputPort], preserveIndexForOuterLookups: Boolean): ErrorOr[Graph]
}

object WomGraphMaker {
  // This apply lets us grab an appropriate WomXMaker[A] out of implicit scope like "val maker = WomXMaker[A]"
  // eg used in the implicit class below.
  def apply[A](implicit maker: WomGraphMaker[A]): WomGraphMaker[A] = maker

  // The restriction [A: WomXMaker] is scala syntax magic for "if there exists in scope a WomXMaker for A"
  implicit class CanMakeGraph[A: WomGraphMaker](val a: A) {
    def toWomGraph(includeGraphNodes: Set[GraphNode],
                   outerLookup: Map[String, GraphNodePort.OutputPort],
                   preserveIndexForOuterLookups: Boolean) =
      WomGraphMaker[A].toWomGraph(a, includeGraphNodes, outerLookup, preserveIndexForOuterLookups)
  }
}
