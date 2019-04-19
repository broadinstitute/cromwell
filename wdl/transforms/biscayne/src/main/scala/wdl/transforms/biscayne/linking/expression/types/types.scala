package wdl.transforms.biscayne.linking.expression

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.TypeEvaluator
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.transforms.base.linking.expression.types.BinaryOperatorEvaluators._
import wdl.transforms.base.linking.expression.types.EngineFunctionEvaluators._
import wdl.transforms.base.linking.expression.types.LiteralEvaluators._
import wdl.transforms.base.linking.expression.types.LookupEvaluators._
import wdl.transforms.base.linking.expression.types.TernaryIfEvaluator._
import wdl.transforms.base.linking.expression.types.UnaryOperatorEvaluators._
import wom.types.WomType
import wdl.transforms.biscayne.linking.expression.types.BiscayneTypeEvaluators._

package object types {
  implicit val expressionTypeEvaluator: TypeEvaluator[ExpressionElement] = new TypeEvaluator[ExpressionElement] {
    override def evaluateType(a: ExpressionElement, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit typeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] = {

      a match {
        // Literals:
        case a: PrimitiveLiteralExpressionElement => a.evaluateType(linkedValues)(typeEvaluator)
        case a: NoneLiteralElement.type => a.evaluateType(linkedValues)(typeEvaluator)

        case a: StringLiteral => a.evaluateType(linkedValues)(typeEvaluator)
        case a: StringExpression => a.evaluateType(linkedValues)(typeEvaluator)
        case a: ObjectLiteral => a.evaluateType(linkedValues)(typeEvaluator)
        case a: MapLiteral => a.evaluateType(linkedValues)(typeEvaluator)
        case a: ArrayLiteral => a.evaluateType(linkedValues)(typeEvaluator)
        case a: PairLiteral => a.evaluateType(linkedValues)(typeEvaluator)

        // Lookups and member accesses:
        case a: IdentifierLookup => a.evaluateType(linkedValues)(typeEvaluator)
        case a: ExpressionMemberAccess => a.evaluateType(linkedValues)(typeEvaluator)
        case a: IdentifierMemberAccess => a.evaluateType(linkedValues)(typeEvaluator)
        case a: IndexAccess => a.evaluateType(linkedValues)(typeEvaluator)

        // Unary operators:
        case a: UnaryNegation => a.evaluateType(linkedValues)(typeEvaluator)
        case a: UnaryPlus => a.evaluateType(linkedValues)(typeEvaluator)
        case a: LogicalNot => a.evaluateType(linkedValues)(typeEvaluator)

        // Binary operators (at some point we might want to split these into separate cases):
        case a: LogicalOr => a.evaluateType(linkedValues)(typeEvaluator)
        case a: LogicalAnd => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Equals => a.evaluateType(linkedValues)(typeEvaluator)
        case a: NotEquals => a.evaluateType(linkedValues)(typeEvaluator)
        case a: LessThan => a.evaluateType(linkedValues)(typeEvaluator)
        case a: LessThanOrEquals => a.evaluateType(linkedValues)(typeEvaluator)
        case a: GreaterThan => a.evaluateType(linkedValues)(typeEvaluator)
        case a: GreaterThanOrEquals => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Add => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Subtract => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Multiply => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Divide => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Remainder => a.evaluateType(linkedValues)(typeEvaluator)

        case a: TernaryIf => a.evaluateType(linkedValues)(typeEvaluator)

        // Engine functions:
        case a: ReadLines => a.evaluateType(linkedValues)(typeEvaluator)
        case a: ReadTsv => a.evaluateType(linkedValues)(typeEvaluator)
        case a: ReadMap => a.evaluateType(linkedValues)(typeEvaluator)
        case a: ReadObject => a.evaluateType(linkedValues)(typeEvaluator)
        case a: ReadObjects => a.evaluateType(linkedValues)(typeEvaluator)
        case a: ReadJson => a.evaluateType(linkedValues)(typeEvaluator)
        case a: ReadInt => a.evaluateType(linkedValues)(typeEvaluator)
        case a: ReadString => a.evaluateType(linkedValues)(typeEvaluator)
        case a: ReadFloat => a.evaluateType(linkedValues)(typeEvaluator)
        case a: ReadBoolean => a.evaluateType(linkedValues)(typeEvaluator)
        case a: WriteLines => a.evaluateType(linkedValues)(typeEvaluator)
        case a: WriteTsv => a.evaluateType(linkedValues)(typeEvaluator)
        case a: WriteMap => a.evaluateType(linkedValues)(typeEvaluator)
        case a: WriteObject => a.evaluateType(linkedValues)(typeEvaluator)
        case a: WriteObjects => a.evaluateType(linkedValues)(typeEvaluator)
        case a: WriteJson => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Range => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Transpose => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Length => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Flatten => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Prefix => a.evaluateType(linkedValues)(typeEvaluator)
        case a: SelectFirst => a.evaluateType(linkedValues)(typeEvaluator)
        case a: SelectAll => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Defined => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Floor => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Ceil => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Round => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Glob => a.evaluateType(linkedValues)(typeEvaluator)

        case a: Size => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Basename => a.evaluateType(linkedValues)(typeEvaluator)

        case a: Zip => a.evaluateType(linkedValues)(typeEvaluator)
        case a: Cross => a.evaluateType(linkedValues)(typeEvaluator)

        case a: Sub => a.evaluateType(linkedValues)(typeEvaluator)

        case a: StdoutElement.type => a.evaluateType(linkedValues)(typeEvaluator)
        case a: StderrElement.type => a.evaluateType(linkedValues)(typeEvaluator)

        case a: Keys => a.evaluateType(linkedValues)(typeEvaluator)
        case a: AsMap => a.evaluateType(linkedValues)(typeEvaluator)
        case a: AsPairs => a.evaluateType(linkedValues)(typeEvaluator)
        case a: CollectByKey => a.evaluateType(linkedValues)(typeEvaluator)

        case other => s"Unable to process ${other.getClass.getSimpleName}: No evaluateType exists for that type.".invalidNel
      }
    }
  }
}
