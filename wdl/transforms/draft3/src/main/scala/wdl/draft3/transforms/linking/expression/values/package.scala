package wdl.draft3.transforms.linking.expression

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.model.draft3.graph.expression.{EvaluatedValue, ForCommandInstantiationOptions, ValueEvaluator}
import wdl.transforms.base.linking.expression.values.BinaryOperatorEvaluators._
import wdl.transforms.base.linking.expression.values.EngineFunctionEvaluators._
import wdl.transforms.base.linking.expression.values.LiteralEvaluators._
import wdl.transforms.base.linking.expression.values.LookupEvaluators._
import wdl.transforms.base.linking.expression.values.TernaryIfEvaluator.ternaryIfEvaluator
import wdl.transforms.base.linking.expression.values.UnaryOperatorEvaluators._
import wom.expression.IoFunctionSet
import wom.values.WomValue
import wdl.transforms.base.wdlom2wdl.WdlWriter.ops._
import wdl.transforms.base.wdlom2wdl.WdlWriterImpl.expressionElementWriter

package object values {

  implicit val expressionEvaluator: ValueEvaluator[ExpressionElement] = new ValueEvaluator[ExpressionElement] {
    override def evaluateValue(a: ExpressionElement,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
    )(implicit valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] =
      a match {
        // Literals:
        case a: PrimitiveLiteralExpressionElement =>
          a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: StringLiteral => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: StringExpression =>
          a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: ObjectLiteral => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: MapLiteral => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: ArrayLiteral => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: PairLiteral => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)

        // Lookups and member accesses:
        case a: IdentifierLookup =>
          a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: ExpressionMemberAccess =>
          a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: IdentifierMemberAccess =>
          a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: IndexAccess => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)

        // Unary operators:
        case a: UnaryNegation => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: UnaryPlus => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: LogicalNot => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)

        // Binary operators (at some point we might want to split these into separate cases):
        case a: LogicalOr => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: LogicalAnd => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Equals => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: NotEquals => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: LessThan => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: LessThanOrEquals =>
          a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: GreaterThan => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: GreaterThanOrEquals =>
          a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Add => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Subtract => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Multiply => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Divide => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Remainder => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)

        case a: TernaryIf => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)

        // Engine functions:
        case a: StdoutElement.type =>
          a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: StderrElement.type =>
          a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)

        case a: ReadLines => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: ReadTsv => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: ReadMap => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: ReadObject => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: ReadObjects => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: ReadJson => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: ReadInt => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: ReadString => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: ReadFloat => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: ReadBoolean => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: WriteLines => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: WriteTsv => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: WriteMap => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: WriteObject => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: WriteObjects => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: WriteJson => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Range => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Transpose => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Length => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Flatten => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Prefix => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: SelectFirst => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: SelectAll => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Defined => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Floor => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Ceil => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Round => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Glob => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)

        case a: Size => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Basename => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)

        case a: Zip => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)
        case a: Cross => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)

        case a: Sub => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)(valueEvaluator)

        case other => s"Unable to process ${other.toWdlV1}: No evaluateValue exists for that type in WDL 1.0".invalidNel
      }
  }
}
