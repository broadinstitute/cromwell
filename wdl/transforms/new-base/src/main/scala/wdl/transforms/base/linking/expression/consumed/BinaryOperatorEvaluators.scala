package wdl.transforms.base.linking.expression.consumed

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.UnlinkedConsumedValueHook
import wdl.model.draft3.graph.ExpressionValueConsumer
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._

object BinaryOperatorEvaluators {
  implicit val logicalOrEvaluator: ExpressionValueConsumer[LogicalOr] = forOperation
  implicit val logicalAndEvaluator: ExpressionValueConsumer[LogicalAnd] = forOperation
  implicit val equalsEvaluator: ExpressionValueConsumer[Equals] = forOperation
  implicit val notEqualsEvaluator: ExpressionValueConsumer[NotEquals] = forOperation
  implicit val lessThanEvaluator: ExpressionValueConsumer[LessThan] = forOperation
  implicit val lessThanOrEqualEvaluator: ExpressionValueConsumer[LessThanOrEquals] = forOperation
  implicit val greaterThanEvaluator: ExpressionValueConsumer[GreaterThan] = forOperation
  implicit val greaterThanOrEqualEvaluator: ExpressionValueConsumer[GreaterThanOrEquals] = forOperation
  implicit val addEvaluator: ExpressionValueConsumer[Add] = forOperation
  implicit val subtractEvaluator: ExpressionValueConsumer[Subtract] = forOperation
  implicit val multiplyEvaluator: ExpressionValueConsumer[Multiply] = forOperation
  implicit val divideEvaluator: ExpressionValueConsumer[Divide] = forOperation
  implicit val remainderEvaluator: ExpressionValueConsumer[Remainder] = forOperation

  private def forOperation[A <: BinaryOperation] = new ExpressionValueConsumer[A] {
    override def expressionConsumedValueHooks(a: A)(implicit
      expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
    ): Set[UnlinkedConsumedValueHook] =
      a.left.expressionConsumedValueHooks ++ a.right.expressionConsumedValueHooks
  }
}
