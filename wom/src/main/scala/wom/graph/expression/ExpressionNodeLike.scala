package wom.graph.expression

import common.Checked
import common.validation.ErrorOr.ErrorOr
import wom.expression.IoFunctionSet
import wom.graph.GraphNode
import wom.graph.GraphNodePort.OutputPort
import wom.values.WomValue

/**
  * Trait for nodes that can be evaluated by the engine
  */
trait ExpressionNodeLike extends GraphNode {
  def evaluate(outputPortLookup: OutputPort => ErrorOr[WomValue],
               ioFunctionSet: IoFunctionSet
  ): Checked[Map[OutputPort, WomValue]]
}
