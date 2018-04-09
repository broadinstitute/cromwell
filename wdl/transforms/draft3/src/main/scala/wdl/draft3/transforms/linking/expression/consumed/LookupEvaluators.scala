package wdl.draft3.transforms.linking.expression.consumed

import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.graph.{ExpressionValueConsumer, UnlinkedCallOutputOrIdentifierAndMemberAccessHook, UnlinkedConsumedValueHook, UnlinkedIdentifierHook}

object LookupEvaluators {

  implicit val identifierLookupUnlinkedValueConsumer: ExpressionValueConsumer[IdentifierLookup] = new ExpressionValueConsumer[IdentifierLookup] {
    override def expressionConsumedValueHooks(a: IdentifierLookup): Set[UnlinkedConsumedValueHook] =
      Set[UnlinkedConsumedValueHook](UnlinkedIdentifierHook(a.identifier))
  }

  implicit val identifierMemberAccessUnlinkedValueConsumer: ExpressionValueConsumer[IdentifierMemberAccess] = new ExpressionValueConsumer[IdentifierMemberAccess] {
    override def expressionConsumedValueHooks(a: IdentifierMemberAccess): Set[UnlinkedConsumedValueHook] =
      Set[UnlinkedConsumedValueHook](UnlinkedCallOutputOrIdentifierAndMemberAccessHook(a.first, a.second))
  }

  implicit val expressionMemberAccessUnlinkedValueConsumer: ExpressionValueConsumer[ExpressionMemberAccess] = new ExpressionValueConsumer[ExpressionMemberAccess] {
    override def expressionConsumedValueHooks(a: ExpressionMemberAccess): Set[UnlinkedConsumedValueHook] = a.expression.expressionConsumedValueHooks
  }
}
