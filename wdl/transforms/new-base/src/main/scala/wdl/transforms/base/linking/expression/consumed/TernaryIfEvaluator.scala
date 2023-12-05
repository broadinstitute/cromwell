package wdl.transforms.base.linking.expression.consumed

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.graph.{ExpressionValueConsumer, UnlinkedConsumedValueHook}

object TernaryIfEvaluator {

  implicit val ternaryIfUnlinkedValueConsumer: ExpressionValueConsumer[TernaryIf] =
    new ExpressionValueConsumer[TernaryIf] {
      override def expressionConsumedValueHooks(
        a: TernaryIf
      )(implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): Set[UnlinkedConsumedValueHook] =
        a.condition.expressionConsumedValueHooks ++ a.ifTrue.expressionConsumedValueHooks ++ a.ifFalse.expressionConsumedValueHooks
    }
}
