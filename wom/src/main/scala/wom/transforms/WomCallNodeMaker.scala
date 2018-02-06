package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.graph.CallNode.CallNodeAndNewNodes
import wom.graph.GraphNodePort

trait WomCallNodeMaker[A] {
  def toWomCallNode(a: A, localLookup: Map[String, GraphNodePort.OutputPort], outerLookup: Map[String, GraphNodePort.OutputPort], preserveIndexForOuterLookups: Boolean): ErrorOr[CallNodeAndNewNodes]
}

object WomCallNodeMaker {
  // This apply lets us grab an appropriate WomXMaker[A] out of implicit scope like "val maker = WomXMaker[A]"
  // eg used in the implicit class below.
  def apply[A](implicit maker: WomCallNodeMaker[A]): WomCallNodeMaker[A] = maker

  // The restriction [A: WomXMaker] is scala syntax magic for "if there exists in scope a WomXMaker for A"
  implicit class CanMakeCallNode[A: WomCallNodeMaker](val a: A) {
    def toWomCallNode(localLookup: Map[String, GraphNodePort.OutputPort],
                      outerLookup: Map[String, GraphNodePort.OutputPort],
                      preserveIndexForOuterLookups: Boolean) =
      WomCallNodeMaker[A].toWomCallNode(a, localLookup, outerLookup, preserveIndexForOuterLookups)
  }
}
