package wdl4s.wom.graph

import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.wom.expression.WomExpression
import wdl4s.wom.graph.GraphNode.GraphNodeSetter
import wdl4s.wom.graph.GraphNodePort.OutputPort

case class GraphNodeInputExpression(inputName: String, expression: WomExpression, inputMapping: Map[String, OutputPort]) {
  private[graph] def instantiateExpression(graphNodeSetter: GraphNodeSetter): ErrorOr[(String, InstantiatedExpression)] = InstantiatedExpression.linkWithInputs(graphNodeSetter, expression, inputMapping) map { (inputName, _) }
}
