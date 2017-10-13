package wom.graph


import lenthall.validation.ErrorOr.ErrorOr
import shapeless.Coproduct
import wom.expression.WomExpression
import wom.graph.CallNode.InputDefinitionPointer
import wom.graph.GraphNodePort.{GraphNodeOutputPort, OutputPort}

final case class ExpressionNode(override val identifier: WomIdentifier, instantiatedExpression: InstantiatedExpression) extends GraphNode {

  val womType = instantiatedExpression.womReturnType
  val singleExpressionOutputPort = GraphNodeOutputPort(localName, womType, this)

  lazy val inputDefinitionPointer = Coproduct[InputDefinitionPointer](singleExpressionOutputPort: OutputPort)

  override val inputPorts = instantiatedExpression.inputPorts
  override val outputPorts: Set[GraphNodePort.OutputPort] = Set(singleExpressionOutputPort)
}

object ExpressionNode {
  def linkWithInputs(nodeIdentifier: WomIdentifier, expression: WomExpression, inputMapping: Map[String, OutputPort]): ErrorOr[ExpressionNode] = {

    println("calling form expression node")
    InstantiatedExpression.instantiateExpressionForNode(ExpressionNode.apply)(nodeIdentifier, expression, inputMapping)
  }
}
