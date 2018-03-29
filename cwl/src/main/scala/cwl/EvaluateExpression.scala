package cwl

import common.validation.ErrorOr.ErrorOr
import cwl.ExpressionEvaluator.{ECMAScriptExpression, ECMAScriptFunction}
import shapeless.Poly1
import wom.values.WomValue

object EvaluateExpression extends Poly1 {
  implicit val script: Case.Aux[ECMAScriptExpression, ParameterContext => ErrorOr[WomValue]] = {
    at {
      ExpressionEvaluator.evalExpression
    }
  }

  implicit val function: Case.Aux[ECMAScriptFunction, ParameterContext => ErrorOr[WomValue]] = {
    at {
      ExpressionEvaluator.evalFunction
    }
  }
}
