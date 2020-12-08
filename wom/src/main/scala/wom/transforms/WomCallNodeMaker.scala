package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.graph.CallNode.CallNodeAndNewNodes
import wom.graph.GraphNodePort
import simulacrum._

@typeclass
trait WomCallNodeMaker[A] {
  def toWomCallNode(a: A, localLookup: Map[String, GraphNodePort.OutputPort], outerLookup: Map[String, GraphNodePort.OutputPort], preserveIndexForOuterLookups: Boolean, inASubworkflow: Boolean): ErrorOr[CallNodeAndNewNodes]
}
