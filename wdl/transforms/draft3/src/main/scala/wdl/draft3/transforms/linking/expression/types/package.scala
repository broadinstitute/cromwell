package wdl.draft3.transforms.linking.expression

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.model.draft3.graph.expression.TypeEvaluator
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wdl.draft3.transforms.linking.expression.types.LookupEvaluators._
import wdl.draft3.transforms.linking.expression.types.LiteralEvaluators._
import wdl.draft3.transforms.linking.expression.types.UnaryOperatorEvaluators._
import wdl.draft3.transforms.linking.expression.types.BinaryOperatorEvaluators._
import wdl.draft3.transforms.linking.expression.types.TernaryIfEvaluator._
import wom.types.WomType

package object types {
  implicit val expressionTypeEvaluator: TypeEvaluator[ExpressionElement] = new TypeEvaluator[ExpressionElement] {
    override def evaluateType(a: ExpressionElement, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {

      a match {
        // Literals:
        case a: PrimitiveLiteralExpressionElement => a.evaluateType(linkedValues)
        case a: StringLiteral => a.evaluateType(linkedValues)
        case a: ObjectLiteral => a.evaluateType(linkedValues)
        case a: MapLiteral => a.evaluateType(linkedValues)
        case a: ArrayLiteral => a.evaluateType(linkedValues)
        case a: PairLiteral => a.evaluateType(linkedValues)

        // Lookups and member accesses:
        case a: IdentifierLookup => a.evaluateType(linkedValues)
        case a: ExpressionMemberAccess => a.evaluateType(linkedValues)
        case a: IdentifierMemberAccess => a.evaluateType(linkedValues)

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

        case other => s"Unable to process ${other.getClass.getSimpleName}: No evaluateValue exists for that type.".invalidNel
      }
    }
  }
}
