package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.graph.CallNode.CallNodeAndNewNodes
import wom.graph.GraphNodePort
import simulacrum._
import scala.language.implicitConversions

@typeclass
trait WomCallNodeMaker[A] {
  @op("toWomCallNode")
  def toWomCallNode(a: A, localLookup: Map[String, GraphNodePort.OutputPort], outerLookup: Map[String, GraphNodePort.OutputPort], preserveIndexForOuterLookups: Boolean): ErrorOr[CallNodeAndNewNodes]
}
