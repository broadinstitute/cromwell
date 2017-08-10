package wdl4s.wom.graph

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import cats.data.Validated.Valid
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.wdl.types.WdlType
import wdl4s.wom.expression.WomExpression
import wdl4s.wom.graph.GraphNodePort.{ConnectedInputPort, InputPort, OutputPort}

class InstantiatedExpression private(val expression: WomExpression, val womReturnType: WdlType, val inputMapping: Map[String, InputPort]) {
  val inputPorts = inputMapping.values.toSet
}

object InstantiatedExpression {

  private[graph] def instantiateExpressionForNode[ExpressionBasedNode <: GraphNode](nodeConstructor: (String, InstantiatedExpression) => ExpressionBasedNode)(name: String, expression: WomExpression, inputMapping: Map[String, OutputPort]): ErrorOr[ExpressionBasedNode] = {
    val graphNodeSetter = new GraphNode.GraphNodeSetter()

    for {
      linkedInputs <- InstantiatedExpression.linkWithInputs(graphNodeSetter, expression, inputMapping)
      expressionNode = nodeConstructor(name, linkedInputs)
      _ = graphNodeSetter._graphNode = expressionNode
    } yield expressionNode
  }

  def linkWithInputs(graphNodeSetter: GraphNode.GraphNodeSetter, expression: WomExpression, inputMapping: Map[String, OutputPort]): ErrorOr[InstantiatedExpression] = {
    def linkInput(input: String): ErrorOr[(String, InputPort)] = if (inputMapping.contains(input)) {
      val upstreamPort = inputMapping(input)
      Valid((input, ConnectedInputPort(input, upstreamPort.womType, upstreamPort, graphNodeSetter.get)))
    } else {
      s"Expression cannot be connected without the input $input.".invalidNel
    }

    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
    for {
      linkedInputList <- expression.inputs.toList traverse linkInput
      linkedInputs = linkedInputList.toMap
      inputTypes = linkedInputs map { case (k, v) => k -> v.womType }
      evaluatedType <- expression.evaluateType(inputTypes)
    } yield new InstantiatedExpression(expression, evaluatedType, linkedInputs)
  }
}