package wdl.draft3.transforms.linking.expression

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.ValueEvaluator
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
    override def evaluateValue(a: ExpressionElement, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {

      a match {
        // Literals:
        case a: PrimitiveLiteralExpressionElement => a.evaluateValue(inputs, ioFunctionSet)
        case a: StringLiteral => a.evaluateValue(inputs, ioFunctionSet)
        case a: ObjectLiteral => a.evaluateValue(inputs, ioFunctionSet)
        case a: MapLiteral => a.evaluateValue(inputs, ioFunctionSet)
        case a: ArrayLiteral => a.evaluateValue(inputs, ioFunctionSet)
        case a: PairLiteral => a.evaluateValue(inputs, ioFunctionSet)

        // Lookups and member accesses:
        case a: IdentifierLookup => a.evaluateValue(inputs, ioFunctionSet)
        case a: ExpressionMemberAccess => a.evaluateValue(inputs, ioFunctionSet)
        case a: IdentifierMemberAccess => a.evaluateValue(inputs, ioFunctionSet)

        // Unary operators:
        case a: UnaryNegation => a.evaluateValue(inputs, ioFunctionSet)
        case a: UnaryPlus => a.evaluateValue(inputs, ioFunctionSet)
        case a: LogicalNot => a.evaluateValue(inputs, ioFunctionSet)

        // Binary operators (at some point we might want to split these into separate cases):
        case a: LogicalOr => a.evaluateValue(inputs, ioFunctionSet)
        case a: LogicalAnd => a.evaluateValue(inputs, ioFunctionSet)
        case a: Equals => a.evaluateValue(inputs, ioFunctionSet)
        case a: NotEquals => a.evaluateValue(inputs, ioFunctionSet)
        case a: LessThan => a.evaluateValue(inputs, ioFunctionSet)
        case a: LessThanOrEquals => a.evaluateValue(inputs, ioFunctionSet)
        case a: GreaterThan => a.evaluateValue(inputs, ioFunctionSet)
        case a: GreaterThanOrEquals => a.evaluateValue(inputs, ioFunctionSet)
        case a: Add => a.evaluateValue(inputs, ioFunctionSet)
        case a: Subtract => a.evaluateValue(inputs, ioFunctionSet)
        case a: Multiply => a.evaluateValue(inputs, ioFunctionSet)
        case a: Divide => a.evaluateValue(inputs, ioFunctionSet)
        case a: Remainder => a.evaluateValue(inputs, ioFunctionSet)

        case a: TernaryIf => a.evaluateValue(inputs, ioFunctionSet)

        // Engine functions:
        case StdoutElement => StdoutElement.evaluateValue(inputs, ioFunctionSet)
        case StderrElement => StderrElement.evaluateValue(inputs, ioFunctionSet)

        case a: ReadLines => a.evaluateValue(inputs, ioFunctionSet)
        case a: ReadTsv => a.evaluateValue(inputs, ioFunctionSet)
        case a: ReadMap => a.evaluateValue(inputs, ioFunctionSet)
        case a: ReadObject => a.evaluateValue(inputs, ioFunctionSet)
        case a: ReadObjects => a.evaluateValue(inputs, ioFunctionSet)
        case a: ReadJson => a.evaluateValue(inputs, ioFunctionSet)
        case a: ReadInt => a.evaluateValue(inputs, ioFunctionSet)
        case a: ReadString => a.evaluateValue(inputs, ioFunctionSet)
        case a: ReadFloat => a.evaluateValue(inputs, ioFunctionSet)
        case a: ReadBoolean => a.evaluateValue(inputs, ioFunctionSet)
        case a: WriteLines => a.evaluateValue(inputs, ioFunctionSet)
        case a: WriteTsv => a.evaluateValue(inputs, ioFunctionSet)
        case a: WriteMap => a.evaluateValue(inputs, ioFunctionSet)
        case a: WriteObject => a.evaluateValue(inputs, ioFunctionSet)
        case a: WriteObjects => a.evaluateValue(inputs, ioFunctionSet)
        case a: WriteJson => a.evaluateValue(inputs, ioFunctionSet)
        case a: Range => a.evaluateValue(inputs, ioFunctionSet)
        case a: Transpose => a.evaluateValue(inputs, ioFunctionSet)
        case a: Length => a.evaluateValue(inputs, ioFunctionSet)
        case a: Flatten => a.evaluateValue(inputs, ioFunctionSet)
        case a: Prefix => a.evaluateValue(inputs, ioFunctionSet)
        case a: SelectFirst => a.evaluateValue(inputs, ioFunctionSet)
        case a: SelectAll => a.evaluateValue(inputs, ioFunctionSet)
        case a: Defined => a.evaluateValue(inputs, ioFunctionSet)
        case a: Floor => a.evaluateValue(inputs, ioFunctionSet)
        case a: Ceil => a.evaluateValue(inputs, ioFunctionSet)
        case a: Round => a.evaluateValue(inputs, ioFunctionSet)

        case a: Size => a.evaluateValue(inputs, ioFunctionSet)
        case a: Basename => a.evaluateValue(inputs, ioFunctionSet)

        case a: Zip => a.evaluateValue(inputs, ioFunctionSet)
        case a: Cross => a.evaluateValue(inputs, ioFunctionSet)

        case a: Sub => a.evaluateValue(inputs, ioFunctionSet)

        case other => s"Unable to process ${other.getClass.getSimpleName}: No evaluateValue exists for that type.".invalidNel
      }
    }
  }
}
