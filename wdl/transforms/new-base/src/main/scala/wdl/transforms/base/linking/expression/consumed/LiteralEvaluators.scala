package wdl.transforms.base.linking.expression.consumed

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.graph.{ExpressionValueConsumer, UnlinkedConsumedValueHook}

object LiteralEvaluators {

  implicit val expressionElementSetUnlinkedValueConsumer: ExpressionValueConsumer[Set[ExpressionElement]] = new ExpressionValueConsumer[Set[ExpressionElement]] {
    override def expressionConsumedValueHooks(elements: Set[ExpressionElement]): Set[UnlinkedConsumedValueHook] =
      elements.flatMap { e: ExpressionElement => e.expressionConsumedValueHooks }
  }
  implicit val objectLiteralUnlinkedValueConsumer: ExpressionValueConsumer[ObjectLiteral] = new ExpressionValueConsumer[ObjectLiteral] {
    override def expressionConsumedValueHooks(o: ObjectLiteral): Set[UnlinkedConsumedValueHook] =
      o.elements.values.toSet[ExpressionElement].expressionConsumedValueHooks
  }

  implicit val mapLiteralUnlinkedValueConsumer: ExpressionValueConsumer[MapLiteral] = new ExpressionValueConsumer[MapLiteral] {
    override def expressionConsumedValueHooks(m: MapLiteral): Set[UnlinkedConsumedValueHook] =
      m.elements.keys.toSet[ExpressionElement].expressionConsumedValueHooks ++
        m.elements.values.toSet[ExpressionElement].expressionConsumedValueHooks
  }

  implicit val pairLiteralUnlinkedValueConsumer: ExpressionValueConsumer[PairLiteral] = new ExpressionValueConsumer[PairLiteral] {
    override def expressionConsumedValueHooks(p: PairLiteral): Set[UnlinkedConsumedValueHook] =
      p.left.expressionConsumedValueHooks ++ p.right.expressionConsumedValueHooks
  }

  implicit val arrayLiteralUnlinkedValueConsumer: ExpressionValueConsumer[ArrayLiteral] = new ExpressionValueConsumer[ArrayLiteral] {
    override def expressionConsumedValueHooks(a: ArrayLiteral): Set[UnlinkedConsumedValueHook] =
      a.elements.toSet[ExpressionElement].expressionConsumedValueHooks
  }

  implicit val stringExpressionUnlinkedValueConsumer: ExpressionValueConsumer[StringExpression] = new ExpressionValueConsumer[StringExpression] {
    override def expressionConsumedValueHooks(a: StringExpression): Set[UnlinkedConsumedValueHook] =
      a.pieces.flatMap {
        case StringLiteral(_) => List.empty
        case StringPlaceholder(expr) => expr.expressionConsumedValueHooks.toList
      }.toSet
  }
}
