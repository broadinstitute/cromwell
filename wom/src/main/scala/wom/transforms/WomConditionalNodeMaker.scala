package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.graph.ConditionalNode.ConditionalNodeWithNewNodes
import wom.graph.GraphNodePort
import wom.graph.GraphNodePort.OutputPort
import simulacrum._

@typeclass
trait WomConditionalNodeMaker[A] {
  def toWomConditionalNode(a: A, localLookup: Map[String, GraphNodePort.OutputPort], outerLookup: Map[String, OutputPort], preserveIndexForOuterLookups: Boolean, inASubworkflow: Boolean): ErrorOr[ConditionalNodeWithNewNodes]
}
