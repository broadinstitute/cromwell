package wdl.draft3.transforms.linking.expression

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.{UnlinkedCallOutputOrIdentifierAndMemberAccessHook, UnlinkedConsumedValueHook, UnlinkedIdentifierHook, ExpressionValueConsumer}
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.elements.ExpressionElement._
import wdl.draft3.transforms.linking.expression.consumed.BinaryOperatorEvaluators._
import wdl.draft3.transforms.linking.expression.consumed.EngineFunctionEvaluators._

package object consumed {

  implicit val expressionElementSetUnlinkedValueConsumer: ExpressionValueConsumer[Set[ExpressionElement]] = new ExpressionValueConsumer[Set[ExpressionElement]] {
    override def expressionConsumedValueHooks(elements: Set[ExpressionElement]): Set[UnlinkedConsumedValueHook] =
      elements.flatMap { e: ExpressionElement => e.expressionConsumedValueHooks }
  }

  implicit val identifierLookupUnlinkedValueConsumer: ExpressionValueConsumer[IdentifierLookup] = new ExpressionValueConsumer[IdentifierLookup] {
    override def expressionConsumedValueHooks(a: IdentifierLookup): Set[UnlinkedConsumedValueHook] =
      Set[UnlinkedConsumedValueHook](UnlinkedIdentifierHook(a.identifier))
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

  implicit val identifierMemberAccessUnlinkedValueConsumer: ExpressionValueConsumer[IdentifierMemberAccess] = new ExpressionValueConsumer[IdentifierMemberAccess] {
    override def expressionConsumedValueHooks(a: IdentifierMemberAccess): Set[UnlinkedConsumedValueHook] =
      Set[UnlinkedConsumedValueHook](UnlinkedCallOutputOrIdentifierAndMemberAccessHook(a.first, a.second))
  }

  implicit val expressionElementUnlinkedValueConsumer: ExpressionValueConsumer[ExpressionElement] = new ExpressionValueConsumer[ExpressionElement] {
    override def expressionConsumedValueHooks(a: ExpressionElement): Set[UnlinkedConsumedValueHook] = a match {
      case _: PrimitiveLiteralExpressionElement | _: StringLiteral => Set.empty[UnlinkedConsumedValueHook]
      case a: ObjectLiteral => a.expressionConsumedValueHooks
      case a: PairLiteral => a.expressionConsumedValueHooks
      case a: ArrayLiteral => a.expressionConsumedValueHooks
      case a: MapLiteral => a.expressionConsumedValueHooks

      case a: IdentifierLookup => a.expressionConsumedValueHooks
      case a: IdentifierMemberAccess => a.expressionConsumedValueHooks

      // Binary operators (at some point we might want to split these into separate cases):
      case a: LogicalOr => a.expressionConsumedValueHooks
      case a: LogicalAnd => a.expressionConsumedValueHooks
      case a: Equals => a.expressionConsumedValueHooks
      case a: NotEquals => a.expressionConsumedValueHooks
      case a: LessThan => a.expressionConsumedValueHooks
      case a: LessThanOrEquals => a.expressionConsumedValueHooks
      case a: GreaterThan => a.expressionConsumedValueHooks
      case a: GreaterThanOrEquals => a.expressionConsumedValueHooks
      case a: Add => a.expressionConsumedValueHooks
      case a: Subtract => a.expressionConsumedValueHooks
      case a: Multiply => a.expressionConsumedValueHooks
      case a: Divide => a.expressionConsumedValueHooks
      case a: Remainder => a.expressionConsumedValueHooks

      // Engine functions:
      case a: Range => a.expressionConsumedValueHooks

      // TODO fill in other expression types
      case other => throw new Exception(s"Cannot generate consumed values for ExpressionElement ${other.getClass.getSimpleName}")
    }
  }
}
