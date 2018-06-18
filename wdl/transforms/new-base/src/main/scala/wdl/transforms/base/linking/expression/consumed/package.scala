package wdl.transforms.base.linking.expression

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.{UnlinkedConsumedValueHook, ExpressionValueConsumer}
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.elements.ExpressionElement._
import wdl.transforms.base.linking.expression.consumed.UnaryOperatorEvaluators._
import wdl.transforms.base.linking.expression.consumed.BinaryOperatorEvaluators._
import wdl.transforms.base.linking.expression.consumed.EngineFunctionEvaluators._
import wdl.transforms.base.linking.expression.consumed.LiteralEvaluators._
import wdl.transforms.base.linking.expression.consumed.LookupEvaluators._
import wdl.transforms.base.linking.expression.consumed.TernaryIfEvaluator._

package object consumed {

  implicit val expressionElementSetUnlinkedValueConsumer: ExpressionValueConsumer[Set[ExpressionElement]] = new ExpressionValueConsumer[Set[ExpressionElement]] {
    override def expressionConsumedValueHooks(elements: Set[ExpressionElement]): Set[UnlinkedConsumedValueHook] =
      elements.flatMap { e: ExpressionElement => e.expressionConsumedValueHooks }
  }

  implicit val kvPairUnlinkedValueConsumer: ExpressionValueConsumer[KvPair] = new ExpressionValueConsumer[KvPair] {
    override def expressionConsumedValueHooks(a: ExpressionElement.KvPair): Set[UnlinkedConsumedValueHook] =
      a.value.expressionConsumedValueHooks
  }

  implicit val expressionElementUnlinkedValueConsumer: ExpressionValueConsumer[ExpressionElement] = new ExpressionValueConsumer[ExpressionElement] {
    override def expressionConsumedValueHooks(a: ExpressionElement): Set[UnlinkedConsumedValueHook] = a match {
      case _: PrimitiveLiteralExpressionElement | _: StringLiteral => Set.empty[UnlinkedConsumedValueHook]
      case a: StringExpression => a.expressionConsumedValueHooks
      case a: ObjectLiteral => a.expressionConsumedValueHooks
      case a: PairLiteral => a.expressionConsumedValueHooks
      case a: ArrayLiteral => a.expressionConsumedValueHooks
      case a: MapLiteral => a.expressionConsumedValueHooks

      // Member access:
      case a: IdentifierLookup => a.expressionConsumedValueHooks
      case a: IdentifierMemberAccess => a.expressionConsumedValueHooks
      case a: ExpressionMemberAccess => a.expressionConsumedValueHooks
      case a: IndexAccess => a.expressionConsumedValueHooks

      // Unary operators:
      case a: UnaryNegation => a.expressionConsumedValueHooks
      case a: UnaryPlus => a.expressionConsumedValueHooks
      case a: LogicalNot => a.expressionConsumedValueHooks

      // Binary operators:
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

      case a: TernaryIf => a.expressionConsumedValueHooks

      // Engine functions:
      case StdoutElement => StdoutElement.expressionConsumedValueHooks
      case StderrElement => StderrElement.expressionConsumedValueHooks

      case a: ReadLines => a.expressionConsumedValueHooks
      case a: ReadTsv => a.expressionConsumedValueHooks
      case a: ReadMap => a.expressionConsumedValueHooks
      case a: ReadObject => a.expressionConsumedValueHooks
      case a: ReadObjects => a.expressionConsumedValueHooks
      case a: ReadJson => a.expressionConsumedValueHooks
      case a: ReadInt => a.expressionConsumedValueHooks
      case a: ReadString => a.expressionConsumedValueHooks
      case a: ReadFloat => a.expressionConsumedValueHooks
      case a: ReadBoolean => a.expressionConsumedValueHooks
      case a: WriteLines => a.expressionConsumedValueHooks
      case a: WriteTsv => a.expressionConsumedValueHooks
      case a: WriteMap => a.expressionConsumedValueHooks
      case a: WriteObject => a.expressionConsumedValueHooks
      case a: WriteObjects => a.expressionConsumedValueHooks
      case a: WriteJson => a.expressionConsumedValueHooks
      case a: Range => a.expressionConsumedValueHooks
      case a: Transpose => a.expressionConsumedValueHooks
      case a: Length => a.expressionConsumedValueHooks
      case a: Flatten => a.expressionConsumedValueHooks
      case a: Prefix => a.expressionConsumedValueHooks
      case a: SelectFirst => a.expressionConsumedValueHooks
      case a: SelectAll => a.expressionConsumedValueHooks
      case a: Defined => a.expressionConsumedValueHooks
      case a: Floor => a.expressionConsumedValueHooks
      case a: Ceil => a.expressionConsumedValueHooks
      case a: Round => a.expressionConsumedValueHooks
      case a: Glob => a.expressionConsumedValueHooks

      case a: Size => a.expressionConsumedValueHooks
      case a: Basename => a.expressionConsumedValueHooks

      case a: Zip => a.expressionConsumedValueHooks
      case a: Cross => a.expressionConsumedValueHooks

      case a: Sub => a.expressionConsumedValueHooks

      // TODO fill in other expression types
      case other => throw new Exception(s"Cannot generate consumed values for ExpressionElement ${other.getClass.getSimpleName}")
    }
  }
}
