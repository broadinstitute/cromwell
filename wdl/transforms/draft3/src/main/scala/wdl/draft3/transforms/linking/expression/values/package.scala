package wdl.draft3.transforms.linking.expression

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.model.draft3.graph.expression.ValueEvaluator
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wom.expression.IoFunctionSet
import wom.values.WomValue
import wdl.draft3.transforms.linking.expression.values.LiteralEvaluators._
import wdl.draft3.transforms.linking.expression.values.UnaryOperatorEvaluators._
import wdl.draft3.transforms.linking.expression.values.BinaryOperatorEvaluators._
import wdl.draft3.transforms.linking.expression.values.LookupEvaluators._
import wdl.draft3.transforms.linking.expression.values.TernaryIfEvaluator.ternaryIfEvaluator

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



        case other => s"Unable to process ${other.getClass.getSimpleName}: No evaluateValue exists for that type.".invalidNel
      }
    }
  }
}
