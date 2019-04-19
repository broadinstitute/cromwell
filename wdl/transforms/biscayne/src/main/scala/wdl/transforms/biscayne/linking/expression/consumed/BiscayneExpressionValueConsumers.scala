package wdl.transforms.biscayne.linking.expression.consumed

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.{ExpressionValueConsumer, UnlinkedConsumedValueHook}

object BiscayneExpressionValueConsumers {
  implicit val keysExpressionValueConsumer: ExpressionValueConsumer[Keys] = new ExpressionValueConsumer[Keys] {
    override def expressionConsumedValueHooks(a: Keys)(implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): Set[UnlinkedConsumedValueHook] = {
      expressionValueConsumer.expressionConsumedValueHooks(a.param)(expressionValueConsumer)
    }
  }

  implicit val asMapExpressionValueConsumer: ExpressionValueConsumer[AsMap] = new ExpressionValueConsumer[AsMap] {
    override def expressionConsumedValueHooks(a: AsMap)(implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): Set[UnlinkedConsumedValueHook] = {
      expressionValueConsumer.expressionConsumedValueHooks(a.param)(expressionValueConsumer)
    }
  }

  implicit val asPairsExpressionValueConsumer: ExpressionValueConsumer[AsPairs] = new ExpressionValueConsumer[AsPairs] {
    override def expressionConsumedValueHooks(a: AsPairs)(implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): Set[UnlinkedConsumedValueHook] = {
      expressionValueConsumer.expressionConsumedValueHooks(a.param)(expressionValueConsumer)
    }
  }

  implicit val collectByKeyExpressionValueConsumer: ExpressionValueConsumer[CollectByKey] = new ExpressionValueConsumer[CollectByKey] {
    override def expressionConsumedValueHooks(a: CollectByKey)(implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): Set[UnlinkedConsumedValueHook] = {
      expressionValueConsumer.expressionConsumedValueHooks(a.param)(expressionValueConsumer)
    }
  }

  implicit val noneLiteralExpressionValueConsumer: ExpressionValueConsumer[NoneLiteralElement.type] = new ExpressionValueConsumer[NoneLiteralElement.type] {
    override def expressionConsumedValueHooks(a: NoneLiteralElement.type)(implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): Set[UnlinkedConsumedValueHook] = {
      // None literals consume no values:
      Set.empty[UnlinkedConsumedValueHook]
    }
  }
}
