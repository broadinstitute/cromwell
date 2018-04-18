package wdl.draft3.transforms.linking.expression.consumed

import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.graph.{ExpressionValueConsumer, UnlinkedCallOutputOrIdentifierAndMemberAccessHook, UnlinkedConsumedValueHook, UnlinkedIdentifierHook}

object TernaryIfEvaluator {

  implicit val ternaryIfUnlinkedValueConsumer: ExpressionValueConsumer[TernaryIf] = new ExpressionValueConsumer[TernaryIf] {
    override def expressionConsumedValueHooks(a: TernaryIf): Set[UnlinkedConsumedValueHook] =
      a.condition.expressionConsumedValueHooks ++ a.ifTrue.expressionConsumedValueHooks ++ a.ifFalse.expressionConsumedValueHooks
  }
}
