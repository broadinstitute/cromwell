package wdl4s.wom.graph


import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.wom.expression.WomExpression
import wdl4s.wom.graph.GraphNodePort.{GraphNodeOutputPort, OutputPort}

case class ExpressionNode private(name: String, instantiatedExpression: InstantiatedExpression) extends GraphNode {

  val womType = instantiatedExpression.womReturnType
  val singleExpressionOutputPort = GraphNodeOutputPort(name, womType, this)

  override val inputPorts = instantiatedExpression.inputPorts
  override val outputPorts: Set[GraphNodePort.OutputPort] = Set(singleExpressionOutputPort)
}

object ExpressionNode {
  def linkWithInputs(name: String, expression: WomExpression, inputMapping: Map[String, OutputPort]): ErrorOr[ExpressionNode] =
    InstantiatedExpression.instantiateExpressionForNode(ExpressionNode.apply)(name, expression, inputMapping)
}
