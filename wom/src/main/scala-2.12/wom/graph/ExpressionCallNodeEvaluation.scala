package wom.graph

import common.Checked
import common.validation.ErrorOr.ErrorOr
import wom.expression.IoFunctionSet
import wom.graph.GraphNodePort.OutputPort
import wom.values.WomValue

object ExpressionCallNodeEvaluation {
  def evaluate(node: ExpressionCallNode, outputPortLookup: OutputPort => ErrorOr[WomValue], ioFunctionSet: IoFunctionSet): Checked[Map[OutputPort, WomValue]] = {
    for {
      womEvaluatedInputs <- CallNode.resolveAndEvaluateInputs(node, ioFunctionSet, outputPortLookup).toEither
      lookup = womEvaluatedInputs.map({ case (inputDefinition, value) => inputDefinition.name -> value })
      evaluated <- node.callable.evaluateFunction(lookup, ioFunctionSet, node.expressionBasedOutputPorts)
    } yield evaluated
  }
}
