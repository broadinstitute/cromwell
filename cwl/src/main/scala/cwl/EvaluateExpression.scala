package cwl

import common.validation.ErrorOr.ErrorOr
import cwl.ExpressionEvaluator.{ECMAScriptExpression, ECMAScriptFunction}
import shapeless.Poly1
import wom.values.WomValue

object EvaluateExpression extends Poly1 {
  implicit def caseECMAScriptExpression: Case.Aux[ECMAScriptExpression, ParameterContext => ErrorOr[WomValue]] = {
    at {
      ExpressionEvaluator.evalExpression
    }
  }

  implicit def caseECMAScriptFunction: Case.Aux[ECMAScriptFunction, ParameterContext => ErrorOr[WomValue]] = {
    at {
      ExpressionEvaluator.evalFunction
    }
  }
}
