package wdl.transforms.base.linking.expression.consumed

import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.graph.{ExpressionValueConsumer, UnlinkedCallOutputOrIdentifierAndMemberAccessHook, UnlinkedConsumedValueHook, UnlinkedIdentifierHook}

object LookupEvaluators {

  implicit val identifierLookupUnlinkedValueConsumer: ExpressionValueConsumer[IdentifierLookup] =
    (a: IdentifierLookup) => Set[UnlinkedConsumedValueHook](UnlinkedIdentifierHook(a.identifier))

  implicit val identifierMemberAccessUnlinkedValueConsumer: ExpressionValueConsumer[IdentifierMemberAccess] =
    (a: IdentifierMemberAccess) => Set[UnlinkedConsumedValueHook](UnlinkedCallOutputOrIdentifierAndMemberAccessHook(a.first, a.second))

  implicit val expressionMemberAccessUnlinkedValueConsumer: ExpressionValueConsumer[ExpressionMemberAccess] =
    (a: ExpressionMemberAccess) => a.expression.expressionConsumedValueHooks

  implicit val indexAccessUnlinkedValueConsumer: ExpressionValueConsumer[IndexAccess] =
    (a: IndexAccess) => a.expressionElement.expressionConsumedValueHooks ++ a.index.expressionConsumedValueHooks
}
