package wom.graph


import lenthall.validation.ErrorOr.ErrorOr
import shapeless.Coproduct
import wom.expression.WomExpression
import wom.graph.CallNode.InputDefinitionPointer
import wom.graph.GraphNodePort.{GraphNodeOutputPort, OutputPort}

final case class ExpressionNode(override val name: String, instantiatedExpression: InstantiatedExpression) extends GraphNode {

  val womType = instantiatedExpression.womReturnType
  val singleExpressionOutputPort = GraphNodeOutputPort(name, womType, this)

  lazy val inputDefinitionPointer = Coproduct[InputDefinitionPointer](singleExpressionOutputPort: OutputPort)

  override val inputPorts = instantiatedExpression.inputPorts
  override val outputPorts: Set[GraphNodePort.OutputPort] = Set(singleExpressionOutputPort)
}

object ExpressionNode {
  def linkWithInputs(name: String, expression: WomExpression, inputMapping: Map[String, OutputPort]): ErrorOr[ExpressionNode] =
    InstantiatedExpression.instantiateExpressionForNode(ExpressionNode.apply)(name, expression, inputMapping)
}
