package wdl.draft3.transforms.linking.expression

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.{EvaluatedValue, ForCommandInstantiationOptions, ValueEvaluator}
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wom.expression.IoFunctionSet
import wom.values.WomValue
import wdl.draft3.transforms.linking.expression.values.LiteralEvaluators._
import wdl.draft3.transforms.linking.expression.values.UnaryOperatorEvaluators._
import wdl.draft3.transforms.linking.expression.values.BinaryOperatorEvaluators._
import wdl.draft3.transforms.linking.expression.values.LookupEvaluators._
import wdl.draft3.transforms.linking.expression.values.TernaryIfEvaluator.ternaryIfEvaluator
import wdl.draft3.transforms.linking.expression.values.EngineFunctionEvaluators._

package object values {

  implicit val expressionEvaluator: ValueEvaluator[ExpressionElement] = new ValueEvaluator[ExpressionElement] {
    override def evaluateValue(a: ExpressionElement,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {

      a match {
        // Literals:
        case a: PrimitiveLiteralExpressionElement => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: StringLiteral => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: StringExpression => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: ObjectLiteral => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: MapLiteral => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: ArrayLiteral => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: PairLiteral => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)

        // Lookups and member accesses:
        case a: IdentifierLookup => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: ExpressionMemberAccess => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: IdentifierMemberAccess => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: IndexAccess => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)

        // Unary operators:
        case a: UnaryNegation => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: UnaryPlus => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: LogicalNot => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)

        // Binary operators (at some point we might want to split these into separate cases):
        case a: LogicalOr => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: LogicalAnd => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Equals => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: NotEquals => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: LessThan => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: LessThanOrEquals => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: GreaterThan => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: GreaterThanOrEquals => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Add => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Subtract => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Multiply => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Divide => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Remainder => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)

        case a: TernaryIf => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)

        // Engine functions:
        case a: StdoutElement.type => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: StderrElement.type => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)

        case a: ReadLines => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: ReadTsv => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: ReadMap => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: ReadObject => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: ReadObjects => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: ReadJson => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: ReadInt => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: ReadString => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: ReadFloat => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: ReadBoolean => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: WriteLines => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: WriteTsv => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: WriteMap => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: WriteObject => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: WriteObjects => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: WriteJson => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Range => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Transpose => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Length => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Flatten => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Prefix => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: SelectFirst => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: SelectAll => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Defined => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Floor => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Ceil => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Round => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Glob => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)

        case a: Size => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Basename => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)

        case a: Zip => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        case a: Cross => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)

        case a: Sub => a.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)

        case other => s"Unable to process ${other.getClass.getSimpleName}: No evaluateValue exists for that type.".invalidNel
      }
    }
  }
}
