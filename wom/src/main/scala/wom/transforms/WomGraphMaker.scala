package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.graph.{Graph, GraphNode, GraphNodePort}
import simulacrum._
import scala.language.implicitConversions

@typeclass
trait WomGraphMaker[A] {
  @op("toWomGraph")
  def toWomGraph(a: A, includeGraphNodes: Set[GraphNode], outerLookup: Map[String, GraphNodePort.OutputPort], preserveIndexForOuterLookups: Boolean): ErrorOr[Graph]
}
