package wom.graph

import lenthall.validation.ErrorOr.ErrorOr
import wom.expression.WomExpression
import wom.graph.GraphNode.GraphNodeSetter
import wom.graph.GraphNodePort.OutputPort
import wom.types.WdlType

/**
  * This is enough information for a GraphNode to build an InstantiatedExpression for an input.
  * Differences:
  * - This one remembers which input the expression is being assigned to.
  * - InstantiatedExpression has created InputPorts for the expression inputs. This one only has references to OutputPorts.
  */
case class GraphNodeInputExpression(inputName: String, expression: WomExpression, inputMapping: Map[String, OutputPort]) {

  /**
    * Instantiate the expression and connect its input ports to the appropriate graphNode.
    */
  private[graph] def instantiateExpression(graphNodeSetter: GraphNodeSetter): ErrorOr[InstantiatedExpression] = InstantiatedExpression.linkWithInputs(graphNodeSetter, expression, inputMapping)

  private[graph] lazy val evaluateType: ErrorOr[WdlType] = expression.evaluateType(inputMapping.map { case (name, port) => (name, port.womType) })
}
