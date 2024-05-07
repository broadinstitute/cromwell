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
    override def evaluateType(a: ExpressionElement,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      typeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      a match {
        // Literals:
        case a: PrimitiveLiteralExpressionElement => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: NoneLiteralElement.type => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)

        case a: StringLiteral => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: StringExpression => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: ObjectLiteral => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: StructLiteral => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: MapLiteral => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: ArrayLiteral => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: PairLiteral => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)

        // Lookups and member accesses:
        case a: IdentifierLookup => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: ExpressionMemberAccess => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: IdentifierMemberAccess => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: IndexAccess => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)

        // Unary operators:
        case a: UnaryNegation => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: UnaryPlus => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: LogicalNot => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)

        // Binary operators (at some point we might want to split these into separate cases):
        case a: LogicalOr => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: LogicalAnd => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Equals => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: NotEquals => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: LessThan => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: LessThanOrEquals => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: GreaterThan => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: GreaterThanOrEquals => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Add => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Subtract => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Multiply => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Divide => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Remainder => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)

        case a: TernaryIf => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)

        // Engine functions:
        case a: ReadLines => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: ReadTsv => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: ReadMap => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: ReadObject => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: ReadObjects => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: ReadJson => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: ReadInt => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: ReadString => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: ReadFloat => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: ReadBoolean => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: WriteLines => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: WriteTsv => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: WriteMap => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: WriteObject => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: WriteObjects => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: WriteJson => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Range => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Transpose => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Length => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Flatten => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Prefix => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Suffix => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: SelectFirst => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: SelectAll => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Defined => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Floor => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Ceil => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Round => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Glob => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Quote => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: SQuote => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Size => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Basename => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)

        case a: Zip => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Cross => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Unzip => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)

        case a: SubPosix => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)

        case a: StdoutElement.type => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: StderrElement.type => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)

        case a: Keys => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: AsMap => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: AsPairs => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: CollectByKey => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Sep => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)

        case a: Min => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)
        case a: Max => a.evaluateType(linkedValues, typeAliases)(typeEvaluator)

        case other =>
          s"Unable to process ${other.getClass.getSimpleName}: No evaluateType exists for that type.".invalidNel
      }
  }
}
