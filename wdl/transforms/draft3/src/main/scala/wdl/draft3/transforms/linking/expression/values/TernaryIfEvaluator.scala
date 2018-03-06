package wdl.draft3.transforms.linking.expression.values

import cats.syntax.validated._

import common.validation.ErrorOr.{ErrorOr, _}
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.ValueEvaluator
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wom.expression.IoFunctionSet
import wom.values.{WomBoolean, WomValue}

object TernaryIfEvaluator {
  implicit val ternaryIfEvaluator: ValueEvaluator[TernaryIf] = new ValueEvaluator[TernaryIf] {
    override def evaluateValue(a: TernaryIf,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {

      a.condition.evaluateValue(inputs, ioFunctionSet) flatMap {
        case WomBoolean(true) => a.ifTrue.evaluateValue(inputs, ioFunctionSet)
        case WomBoolean(false) => a.ifFalse.evaluateValue(inputs, ioFunctionSet)
        case other => s"Condition should have evaluated to a Boolean but instead got ${other.womType.toDisplayString}".invalidNel
      }
    }
  }
}
