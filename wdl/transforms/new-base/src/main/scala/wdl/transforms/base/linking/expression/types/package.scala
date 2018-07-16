package wdl.transforms.base.linking.expression

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.model.draft3.graph.expression.TypeEvaluator
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wdl.transforms.base.linking.expression.types.LookupEvaluators._
import wdl.transforms.base.linking.expression.types.LiteralEvaluators._
import wdl.transforms.base.linking.expression.types.UnaryOperatorEvaluators._
import wdl.transforms.base.linking.expression.types.BinaryOperatorEvaluators._
import wdl.transforms.base.linking.expression.types.TernaryIfEvaluator._
import wdl.transforms.base.linking.expression.types.EngineFunctionEvaluators._
import wom.types.WomType

package object types {
  implicit val expressionTypeEvaluator: TypeEvaluator[ExpressionElement] = new TypeEvaluator[ExpressionElement] {
    override def evaluateType(a: ExpressionElement, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {

      a match {
        // Literals:
        case a: PrimitiveLiteralExpressionElement => a.evaluateType(linkedValues)
        case a: StringLiteral => a.evaluateType(linkedValues)
        case a: StringExpression => a.evaluateType(linkedValues)
        case a: ObjectLiteral => a.evaluateType(linkedValues)
        case a: MapLiteral => a.evaluateType(linkedValues)
        case a: ArrayLiteral => a.evaluateType(linkedValues)
        case a: PairLiteral => a.evaluateType(linkedValues)

        // Lookups and member accesses:
        case a: IdentifierLookup => a.evaluateType(linkedValues)
        case a: ExpressionMemberAccess => a.evaluateType(linkedValues)
        case a: IdentifierMemberAccess => a.evaluateType(linkedValues)
        case a: IndexAccess => a.evaluateType(linkedValues)

        // Unary operators:
        case a: UnaryNegation => a.evaluateType(linkedValues)
        case a: UnaryPlus => a.evaluateType(linkedValues)
        case a: LogicalNot => a.evaluateType(linkedValues)

        // Binary operators (at some point we might want to split these into separate cases):
        case a: LogicalOr => a.evaluateType(linkedValues)
        case a: LogicalAnd => a.evaluateType(linkedValues)
        case a: Equals => a.evaluateType(linkedValues)
        case a: NotEquals => a.evaluateType(linkedValues)
        case a: LessThan => a.evaluateType(linkedValues)
        case a: LessThanOrEquals => a.evaluateType(linkedValues)
        case a: GreaterThan => a.evaluateType(linkedValues)
        case a: GreaterThanOrEquals => a.evaluateType(linkedValues)
        case a: Add => a.evaluateType(linkedValues)
        case a: Subtract => a.evaluateType(linkedValues)
        case a: Multiply => a.evaluateType(linkedValues)
        case a: Divide => a.evaluateType(linkedValues)
        case a: Remainder => a.evaluateType(linkedValues)

        case a: TernaryIf => a.evaluateType(linkedValues)

        // Engine functions:
        case a: ReadLines => a.evaluateType(linkedValues)
        case a: ReadTsv => a.evaluateType(linkedValues)
        case a: ReadMap => a.evaluateType(linkedValues)
        case a: ReadObject => a.evaluateType(linkedValues)
        case a: ReadObjects => a.evaluateType(linkedValues)
        case a: ReadJson => a.evaluateType(linkedValues)
        case a: ReadInt => a.evaluateType(linkedValues)
        case a: ReadString => a.evaluateType(linkedValues)
        case a: ReadFloat => a.evaluateType(linkedValues)
        case a: ReadBoolean => a.evaluateType(linkedValues)
        case a: WriteLines => a.evaluateType(linkedValues)
        case a: WriteTsv => a.evaluateType(linkedValues)
        case a: WriteMap => a.evaluateType(linkedValues)
        case a: WriteObject => a.evaluateType(linkedValues)
        case a: WriteObjects => a.evaluateType(linkedValues)
        case a: WriteJson => a.evaluateType(linkedValues)
        case a: Range => a.evaluateType(linkedValues)
        case a: Transpose => a.evaluateType(linkedValues)
        case a: Length => a.evaluateType(linkedValues)
        case a: Flatten => a.evaluateType(linkedValues)
        case a: Prefix => a.evaluateType(linkedValues)
        case a: SelectFirst => a.evaluateType(linkedValues)
        case a: SelectAll => a.evaluateType(linkedValues)
        case a: Defined => a.evaluateType(linkedValues)
        case a: Floor => a.evaluateType(linkedValues)
        case a: Ceil => a.evaluateType(linkedValues)
        case a: Round => a.evaluateType(linkedValues)
        case a: Glob => a.evaluateType(linkedValues)

        case a: Size => a.evaluateType(linkedValues)
        case a: Basename => a.evaluateType(linkedValues)

        case a: Zip => a.evaluateType(linkedValues)
        case a: Cross => a.evaluateType(linkedValues)

        case a: Sub => a.evaluateType(linkedValues)

        case other => s"Unable to process ${other.getClass.getSimpleName}: No evaluateType exists for that type.".invalidNel
      }
    }
  }
}
