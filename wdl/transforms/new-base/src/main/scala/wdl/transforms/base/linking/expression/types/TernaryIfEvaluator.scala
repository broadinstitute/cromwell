package wdl.transforms.base.linking.expression.types

import cats.syntax.apply._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.model.draft3.graph.expression.TypeEvaluator
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wom.types.{WomBooleanType, WomType}

object TernaryIfEvaluator {
  implicit val ternaryIfEvaluator: TypeEvaluator[TernaryIf] = new TypeEvaluator[TernaryIf] {
    override def evaluateType(a: TernaryIf, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] = {
      a.condition.evaluateType(linkedValues) flatMap {
        case WomBooleanType =>
          (a.ifTrue.evaluateType(linkedValues): ErrorOr[WomType],
            a.ifFalse.evaluateType(linkedValues): ErrorOr[WomType]) mapN { (tType, fType) => WomType.homogeneousTypeFromTypes(Seq(tType, fType)) }
        case other => s"Condition should have evaluated to a Boolean but instead got ${other.stableName}".invalidNel
      }
    }
  }
}
