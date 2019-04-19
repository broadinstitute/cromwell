package wdl.transforms.biscayne.linking.expression

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.graph.{ExpressionValueConsumer, UnlinkedConsumedValueHook}
import wdl.transforms.base.linking.expression.consumed.BinaryOperatorEvaluators._
import wdl.transforms.base.linking.expression.consumed.EngineFunctionEvaluators._
import wdl.transforms.base.linking.expression.consumed.LiteralEvaluators._
import wdl.transforms.base.linking.expression.consumed.LookupEvaluators._
import wdl.transforms.base.linking.expression.consumed.TernaryIfEvaluator._
import wdl.transforms.base.linking.expression.consumed.UnaryOperatorEvaluators._
import wdl.transforms.biscayne.linking.expression.consumed.BiscayneExpressionValueConsumers._

package object consumed {

  implicit val expressionElementUnlinkedValueConsumer: ExpressionValueConsumer[ExpressionElement] = new ExpressionValueConsumer[ExpressionElement] {
    override def expressionConsumedValueHooks(a: ExpressionElement)(implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): Set[UnlinkedConsumedValueHook] = a match {
      case _: PrimitiveLiteralExpressionElement | _: StringLiteral => Set.empty[UnlinkedConsumedValueHook]
      case a: NoneLiteralElement.type => a.expressionConsumedValueHooks(expressionValueConsumer)

      case a: StringExpression => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: ObjectLiteral => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: PairLiteral => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: ArrayLiteral => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: MapLiteral => a.expressionConsumedValueHooks(expressionValueConsumer)

      // Member access:
      case a: IdentifierLookup => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: IdentifierMemberAccess => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: ExpressionMemberAccess => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: IndexAccess => a.expressionConsumedValueHooks(expressionValueConsumer)

      // Unary operators:
      case a: UnaryNegation => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: UnaryPlus => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: LogicalNot => a.expressionConsumedValueHooks(expressionValueConsumer)

      // Binary operators:
      case a: LogicalOr => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: LogicalAnd => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Equals => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: NotEquals => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: LessThan => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: LessThanOrEquals => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: GreaterThan => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: GreaterThanOrEquals => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Add => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Subtract => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Multiply => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Divide => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Remainder => a.expressionConsumedValueHooks(expressionValueConsumer)

      case a: TernaryIf => a.expressionConsumedValueHooks(expressionValueConsumer)

      // Engine functions:
      case a: StdoutElement.type => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: StderrElement.type => a.expressionConsumedValueHooks(expressionValueConsumer)

      case a: ReadLines => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: ReadTsv => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: ReadMap => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: ReadObject => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: ReadObjects => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: ReadJson => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: ReadInt => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: ReadString => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: ReadFloat => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: ReadBoolean => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: WriteLines => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: WriteTsv => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: WriteMap => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: WriteObject => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: WriteObjects => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: WriteJson => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Range => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Transpose => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Length => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Flatten => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Prefix => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: SelectFirst => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: SelectAll => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Defined => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Floor => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Ceil => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Round => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Glob => a.expressionConsumedValueHooks(expressionValueConsumer)

      case a: Size => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Basename => a.expressionConsumedValueHooks(expressionValueConsumer)

      case a: Zip => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: Cross => a.expressionConsumedValueHooks(expressionValueConsumer)

      case a: Sub => a.expressionConsumedValueHooks(expressionValueConsumer)

      // New WDL biscayne expressions:
      case a: Keys => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: AsMap => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: AsPairs => a.expressionConsumedValueHooks(expressionValueConsumer)
      case a: CollectByKey => a.expressionConsumedValueHooks(expressionValueConsumer)

      case other => throw new Exception(s"Cannot generate consumed values for ExpressionElement ${other.getClass.getSimpleName}")
    }
  }
}
