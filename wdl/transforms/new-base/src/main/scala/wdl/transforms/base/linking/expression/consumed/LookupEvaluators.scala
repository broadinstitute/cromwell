package wdl.transforms.base.linking.expression.consumed

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.graph.{
  ExpressionValueConsumer,
  UnlinkedCallOutputOrIdentifierAndMemberAccessHook,
  UnlinkedConsumedValueHook,
  UnlinkedIdentifierHook
}

object LookupEvaluators {

  implicit val identifierLookupUnlinkedValueConsumer: ExpressionValueConsumer[IdentifierLookup] =
    new ExpressionValueConsumer[IdentifierLookup] {
      override def expressionConsumedValueHooks(a: IdentifierLookup)(implicit
        expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
      ): Set[UnlinkedConsumedValueHook] =
        Set[UnlinkedConsumedValueHook](UnlinkedIdentifierHook(a.identifier))
    }

  implicit val identifierMemberAccessUnlinkedValueConsumer = new ExpressionValueConsumer[IdentifierMemberAccess] {
    override def expressionConsumedValueHooks(a: IdentifierMemberAccess)(implicit
      expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
    ): Set[UnlinkedConsumedValueHook] =
      Set[UnlinkedConsumedValueHook](UnlinkedCallOutputOrIdentifierAndMemberAccessHook(a.first, a.second))
  }

  implicit val expressionMemberAccessUnlinkedValueConsumer = new ExpressionValueConsumer[ExpressionMemberAccess] {
    override def expressionConsumedValueHooks(a: ExpressionMemberAccess)(implicit
      expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
    ): Set[UnlinkedConsumedValueHook] =
      a.expression.expressionConsumedValueHooks
  }

  implicit val indexAccessUnlinkedValueConsumer = new ExpressionValueConsumer[IndexAccess] {
    override def expressionConsumedValueHooks(a: IndexAccess)(implicit
      expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
    ): Set[UnlinkedConsumedValueHook] =
      a.expressionElement.expressionConsumedValueHooks ++ a.index.expressionConsumedValueHooks
  }
}
