package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.graph.GraphNodePort
import wom.graph.GraphNodePort.OutputPort
import wom.graph.ScatterNode.ScatterNodeWithNewNodes

trait WomScatterNodeMaker[A] {
  def toWomScatterNode(a: A, localLookup: Map[String, GraphNodePort.OutputPort], outerLookup: Map[String, OutputPort], preserveIndexForOuterLookups: Boolean): ErrorOr[ScatterNodeWithNewNodes]
}


object WomScatterNodeMaker {
  // This apply lets us grab an appropriate WomXMaker[A] out of implicit scope like "val maker = WomXMaker[A]"
  // eg used in the implicit class below.
  def apply[A](implicit maker: WomScatterNodeMaker[A]): WomScatterNodeMaker[A] = maker

  // The restriction [A: WomXMaker] is scala syntax magic for "if there exists in scope a WomXMaker for A"
  implicit class CanMakeScatterNode[A: WomScatterNodeMaker](val a: A) {
    def toWomScatterNode(localLookup: Map[String, GraphNodePort.OutputPort],
                         outerLookup: Map[String, GraphNodePort.OutputPort],
                         preserveIndexForOuterLookups: Boolean) =
      WomScatterNodeMaker[A].toWomScatterNode(a, localLookup, outerLookup, preserveIndexForOuterLookups)
  }
}

