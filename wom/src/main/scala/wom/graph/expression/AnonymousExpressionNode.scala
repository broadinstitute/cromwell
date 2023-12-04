package wom.graph.expression

import common.validation.ErrorOr.ErrorOr
import wom.expression.WomExpression
import wom.graph.GraphNode.GraphNodeSetter
import wom.graph.GraphNodePort.{InputPort, OutputPort}
import wom.graph.{CommandCallNode, WomIdentifier}
import wom.types.WomType

object AnonymousExpressionNode {
  type AnonymousExpressionConstructor[T] = (WomIdentifier, WomExpression, WomType, Map[String, InputPort]) => T

  def fromInputMapping[T <: AnonymousExpressionNode](identifier: WomIdentifier,
                                                     expression: WomExpression,
                                                     inputMapping: Map[String, OutputPort],
                                                     constructor: AnonymousExpressionConstructor[T]
  ): ErrorOr[T] =
    ExpressionNode.buildFromConstructor(constructor)(identifier, expression, inputMapping)
}

/**
  * An expression node that is purely an internal expression and shouldn't be visible outside the graph
  */
trait AnonymousExpressionNode extends ExpressionNode

case class PlainAnonymousExpressionNode(override val identifier: WomIdentifier,
                                        override val womExpression: WomExpression,
                                        override val womType: WomType,
                                        override val inputMapping: Map[String, InputPort]
) extends ExpressionNode(identifier, womExpression, womType, inputMapping)
    with AnonymousExpressionNode

case class TaskCallInputExpressionNode(override val identifier: WomIdentifier,
                                       override val womExpression: WomExpression,
                                       override val womType: WomType,
                                       override val inputMapping: Map[String, InputPort]
) extends ExpressionNode(identifier, womExpression, womType, inputMapping)
    with AnonymousExpressionNode {

  /**
   * The `GraphNodeSetter` that will have a `TaskCallNode` set into it. This is needed in the `WorkflowExecutionActor`
   * to be able to look up backend mapping for the target call in order to have the correct `IoFunctionSet` to
   * evaluate task call input expressions.
   */
  val taskCallNodeReceivingInput: GraphNodeSetter[CommandCallNode] = new GraphNodeSetter[CommandCallNode]
}
