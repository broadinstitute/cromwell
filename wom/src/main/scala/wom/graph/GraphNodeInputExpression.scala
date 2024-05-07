package wom.graph

import common.validation.ErrorOr.ErrorOr
import wom.expression.WomExpression
import wom.graph.GraphNodePort.OutputPort
import wom.types.WomType

/**
  * This is enough information for a GraphNode to build an InstantiatedExpression for an input.
  * Differences:
  * - This one remembers which input the expression is being assigned to.
  * - InstantiatedExpression has created InputPorts for the expression inputs. This one only has references to OutputPorts.
  */
case class GraphNodeInputExpression(inputName: String,
                                    expression: WomExpression,
                                    inputMapping: Map[String, OutputPort]
) {

  private[graph] lazy val evaluateType: ErrorOr[WomType] = expression.evaluateType(inputMapping.map {
    case (name, port) => (name, port.womType)
  })
}
