package wom.graph

import lenthall.validation.ErrorOr.ErrorOr
import wdl.types.WdlType
import wom.expression.WomExpression
import wom.graph.GraphNodePort.{InputPort, OutputPort}

class InstantiatedExpression (val expression: WomExpression, val womReturnType: WdlType, val inputMapping: Map[String, InputPort]) {
  val inputPorts = inputMapping.values.toSet
}

object InstantiatedExpression {

  private[graph] def instantiateExpressionForNode[ExpressionBasedNode <: GraphNode, Identifier <: WomIdentifier](nodeConstructor: (Identifier, InstantiatedExpression) => ExpressionBasedNode)(nodeIdentifier: Identifier, expression: WomExpression, inputMapping: Map[String, OutputPort]): ErrorOr[ExpressionBasedNode] = {
    val graphNodeSetter = new GraphNode.GraphNodeSetter()

    for {
      linkedInputs <- expression.linkWithInputs(graphNodeSetter, inputMapping)
      expressionNode = nodeConstructor(nodeIdentifier, linkedInputs)
      _ = graphNodeSetter._graphNode = expressionNode
    } yield expressionNode
  }

}