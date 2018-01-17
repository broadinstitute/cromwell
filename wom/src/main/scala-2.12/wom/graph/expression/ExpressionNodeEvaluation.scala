package wom.graph.expression

import common.Checked
import common.validation.ErrorOr.ErrorOr
import wom.expression.IoFunctionSet
import wom.graph.GraphNodePort.OutputPort
import wom.values.WomValue

object ExpressionNodeEvaluation {
  def evaluate(node: ExpressionNode, outputPortLookup: OutputPort => ErrorOr[WomValue], ioFunctionSet: IoFunctionSet): Checked[Map[OutputPort, WomValue]] = {
    import common.validation.ErrorOr._
    for {
      inputs <- node.inputMapping.traverseValues(inputPort => outputPortLookup(inputPort.upstream)).toEither
      evaluated <- node.evaluateAndCoerce(inputs, ioFunctionSet)
    } yield Map(node.singleExpressionOutputPort -> evaluated)
  }
}
