package wdl.transforms.base.linking.expression.consumed

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.graph.{ExpressionValueConsumer, UnlinkedConsumedValueHook}

object UnaryOperatorEvaluators {
  implicit val logicalNotEvaluator: ExpressionValueConsumer[LogicalNot] = forOperation
  implicit val unaryNegationEvaluator: ExpressionValueConsumer[UnaryNegation] = forOperation
  implicit val unaryPlusEvaluator: ExpressionValueConsumer[UnaryPlus] = forOperation

  private def forOperation[A <: UnaryOperation] = new ExpressionValueConsumer[A] {
    override def expressionConsumedValueHooks(a: A)(implicit
      expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
    ): Set[UnlinkedConsumedValueHook] = a.argument.expressionConsumedValueHooks
  }
}
