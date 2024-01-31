package wdl.transforms.cascades.linking.expression.consumed

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.{ExpressionValueConsumer, UnlinkedConsumedValueHook}

object cascadesExpressionValueConsumers {
  implicit val keysExpressionValueConsumer: ExpressionValueConsumer[Keys] = new ExpressionValueConsumer[Keys] {
    override def expressionConsumedValueHooks(a: Keys)(implicit
      expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
    ): Set[UnlinkedConsumedValueHook] =
      expressionValueConsumer.expressionConsumedValueHooks(a.param)(expressionValueConsumer)
  }

  implicit val asMapExpressionValueConsumer: ExpressionValueConsumer[AsMap] = new ExpressionValueConsumer[AsMap] {
    override def expressionConsumedValueHooks(a: AsMap)(implicit
      expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
    ): Set[UnlinkedConsumedValueHook] =
      expressionValueConsumer.expressionConsumedValueHooks(a.param)(expressionValueConsumer)
  }

  implicit val asPairsExpressionValueConsumer: ExpressionValueConsumer[AsPairs] = new ExpressionValueConsumer[AsPairs] {
    override def expressionConsumedValueHooks(a: AsPairs)(implicit
      expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
    ): Set[UnlinkedConsumedValueHook] =
      expressionValueConsumer.expressionConsumedValueHooks(a.param)(expressionValueConsumer)
  }

  implicit val collectByKeyExpressionValueConsumer: ExpressionValueConsumer[CollectByKey] =
    new ExpressionValueConsumer[CollectByKey] {
      override def expressionConsumedValueHooks(a: CollectByKey)(implicit
        expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
      ): Set[UnlinkedConsumedValueHook] =
        expressionValueConsumer.expressionConsumedValueHooks(a.param)(expressionValueConsumer)
    }

  implicit val minExpressionValueConsumer: ExpressionValueConsumer[Min] = new ExpressionValueConsumer[Min] {
    override def expressionConsumedValueHooks(
      a: Min
    )(implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): Set[UnlinkedConsumedValueHook] =
      expressionValueConsumer.expressionConsumedValueHooks(a.arg1)(expressionValueConsumer) ++ expressionValueConsumer
        .expressionConsumedValueHooks(a.arg2)(expressionValueConsumer)
  }

  implicit val maxExpressionValueConsumer: ExpressionValueConsumer[Max] = new ExpressionValueConsumer[Max] {
    override def expressionConsumedValueHooks(
      a: Max
    )(implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): Set[UnlinkedConsumedValueHook] =
      expressionValueConsumer.expressionConsumedValueHooks(a.arg1)(expressionValueConsumer) ++ expressionValueConsumer
        .expressionConsumedValueHooks(a.arg2)(expressionValueConsumer)
  }

  implicit val sepExpressionValueConsumer: ExpressionValueConsumer[Sep] = new ExpressionValueConsumer[Sep] {
    override def expressionConsumedValueHooks(a: Sep)(implicit
      expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
    ): Set[UnlinkedConsumedValueHook] =
      expressionValueConsumer.expressionConsumedValueHooks(a.arg1)(expressionValueConsumer) ++
        expressionValueConsumer.expressionConsumedValueHooks(a.arg2)(expressionValueConsumer)
  }

  implicit val suffixExpressionValueConsumer: ExpressionValueConsumer[Suffix] = new ExpressionValueConsumer[Suffix] {
    override def expressionConsumedValueHooks(a: Suffix)(implicit
      expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
    ): Set[UnlinkedConsumedValueHook] =
      expressionValueConsumer.expressionConsumedValueHooks(a.arg1)(expressionValueConsumer) ++
        expressionValueConsumer.expressionConsumedValueHooks(a.arg2)(expressionValueConsumer)
  }

  implicit val noneLiteralExpressionValueConsumer: ExpressionValueConsumer[NoneLiteralElement.type] =
    new ExpressionValueConsumer[NoneLiteralElement.type] {
      override def expressionConsumedValueHooks(a: NoneLiteralElement.type)(implicit
        expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
      ): Set[UnlinkedConsumedValueHook] =
        // None literals consume no values:
        Set.empty[UnlinkedConsumedValueHook]
    }
}
