package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.graph.{Graph, GraphNode, GraphNodePort}
import simulacrum._

@typeclass
trait WomGraphMaker[A] {
  def toWomGraph(a: A, includeGraphNodes: Set[GraphNode], outerLookup: Map[String, GraphNodePort.OutputPort], preserveIndexForOuterLookups: Boolean, inASubworkflow: Boolean): ErrorOr[Graph]
}
