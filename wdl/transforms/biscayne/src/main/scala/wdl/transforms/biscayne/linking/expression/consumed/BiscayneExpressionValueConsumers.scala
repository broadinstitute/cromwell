package wdl.transforms.biscayne.linking.expression.consumed

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.{ExpressionValueConsumer, UnlinkedConsumedValueHook}

object BiscayneExpressionValueConsumers {
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

  implicit val subPosixExpressionValueConsumer: ExpressionValueConsumer[SubPosix] =
    new ExpressionValueConsumer[SubPosix] {
      override def expressionConsumedValueHooks(a: SubPosix)(implicit
        expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
      ): Set[UnlinkedConsumedValueHook] =
        expressionValueConsumer.expressionConsumedValueHooks(a.arg1)(expressionValueConsumer) ++
          expressionValueConsumer.expressionConsumedValueHooks(a.arg2)(expressionValueConsumer) ++
          expressionValueConsumer.expressionConsumedValueHooks(a.arg3)(expressionValueConsumer)
    }

  implicit val suffixExpressionValueConsumer: ExpressionValueConsumer[Suffix] = new ExpressionValueConsumer[Suffix] {
    override def expressionConsumedValueHooks(a: Suffix)(implicit
      expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
    ): Set[UnlinkedConsumedValueHook] =
      expressionValueConsumer.expressionConsumedValueHooks(a.arg1)(expressionValueConsumer) ++
        expressionValueConsumer.expressionConsumedValueHooks(a.arg2)(expressionValueConsumer)
  }

  implicit val quoteExpressionValueConsumer: ExpressionValueConsumer[Quote] = new ExpressionValueConsumer[Quote] {
    override def expressionConsumedValueHooks(a: Quote)(implicit
      expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
    ): Set[UnlinkedConsumedValueHook] =
      expressionValueConsumer.expressionConsumedValueHooks(a.param)(expressionValueConsumer)
  }

  implicit val sQuoteExpressionValueConsumer: ExpressionValueConsumer[SQuote] = new ExpressionValueConsumer[SQuote] {
    override def expressionConsumedValueHooks(a: SQuote)(implicit
      expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
    ): Set[UnlinkedConsumedValueHook] =
      expressionValueConsumer.expressionConsumedValueHooks(a.param)(expressionValueConsumer)
  }

  implicit val unzipExpressionValueConsumer: ExpressionValueConsumer[Unzip] = new ExpressionValueConsumer[Unzip] {
    override def expressionConsumedValueHooks(a: Unzip)(implicit
      expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
    ): Set[UnlinkedConsumedValueHook] =
      expressionValueConsumer.expressionConsumedValueHooks(a.param)(expressionValueConsumer)
  }

  implicit val noneLiteralExpressionValueConsumer: ExpressionValueConsumer[NoneLiteralElement.type] =
    new ExpressionValueConsumer[NoneLiteralElement.type] {
      override def expressionConsumedValueHooks(a: NoneLiteralElement.type)(implicit
        expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
      ): Set[UnlinkedConsumedValueHook] =
        // None literals consume no values:
        Set.empty[UnlinkedConsumedValueHook]
    }

  implicit val structLiteralExpressionValueConsumer: ExpressionValueConsumer[StructLiteral] =
    new ExpressionValueConsumer[StructLiteral] {
      override def expressionConsumedValueHooks(a: StructLiteral)(implicit
        expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]
      ): Set[UnlinkedConsumedValueHook] =
        a.elements.values
          .flatMap(element => expressionValueConsumer.expressionConsumedValueHooks(element)(expressionValueConsumer))
          .toSet
    }
}
