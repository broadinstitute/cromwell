package wdl.transforms.base.linking.expression.values

import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, _}
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.{EvaluatedValue, ForCommandInstantiationOptions, ValueEvaluator}
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wom.expression.IoFunctionSet
import wom.values.{WomBoolean, WomValue}

object TernaryIfEvaluator {
  implicit val ternaryIfEvaluator: ValueEvaluator[TernaryIf] = new ValueEvaluator[TernaryIf] {
    override def evaluateValue(a: TernaryIf,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {

      a.condition.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions) flatMap {
        case EvaluatedValue(WomBoolean(true), conditionSideEffectFiles) =>
          a.ifTrue.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions).map(result => result.copy(sideEffectFiles = result.sideEffectFiles ++ conditionSideEffectFiles))
        case EvaluatedValue(WomBoolean(false), conditionSideEffectFiles) =>
          a.ifFalse.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions).map(result => result.copy(sideEffectFiles = result.sideEffectFiles ++ conditionSideEffectFiles))
        case other => s"Condition should have evaluated to a Boolean but instead got ${other.value.womType.stableName}".invalidNel
      }
    }
  }
}
