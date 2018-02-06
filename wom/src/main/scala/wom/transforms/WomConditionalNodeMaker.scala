package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.graph.ConditionalNode.ConditionalNodeWithNewNodes
import wom.graph.GraphNodePort
import wom.graph.GraphNodePort.OutputPort

trait WomConditionalNodeMaker[A] {
  def toWomConditionalNode(a: A, localLookup: Map[String, GraphNodePort.OutputPort], outerLookup: Map[String, OutputPort], preserveIndexForOuterLookups: Boolean): ErrorOr[ConditionalNodeWithNewNodes]
}

object WomConditionalNodeMaker {
  // This apply lets us grab an appropriate WomXMaker[A] out of implicit scope like "val maker = WomXMaker[A]"
  // eg used in the implicit class below.
  def apply[A](implicit maker: WomConditionalNodeMaker[A]): WomConditionalNodeMaker[A] = maker

  // The restriction [A: WomXMaker] is scala syntax magic for "if there exists in scope a WomXMaker for A"
  implicit class CanMakeConditionalNode[A: WomConditionalNodeMaker](val a: A) {
    def toWomConditionalNode(localLookup: Map[String, GraphNodePort.OutputPort],
                             outerLookup: Map[String, GraphNodePort.OutputPort],
                             preserveIndexForOuterLookups: Boolean) =
      WomConditionalNodeMaker[A].toWomConditionalNode(a, localLookup, outerLookup, preserveIndexForOuterLookups)
  }
}
