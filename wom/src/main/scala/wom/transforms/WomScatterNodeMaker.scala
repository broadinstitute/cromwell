package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.graph.GraphNodePort
import wom.graph.GraphNodePort.OutputPort
import wom.graph.ScatterNode.ScatterNodeWithNewNodes
import simulacrum._

@typeclass
trait WomScatterNodeMaker[A] {
  def toWomScatterNode(a: A, localLookup: Map[String, GraphNodePort.OutputPort], outerLookup: Map[String, OutputPort], preserveIndexForOuterLookups: Boolean, inASubworkflow: Boolean): ErrorOr[ScatterNodeWithNewNodes]
}
