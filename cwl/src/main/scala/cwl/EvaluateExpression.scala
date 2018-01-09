package cwl

import common.validation.ErrorOr.ErrorOr
import cwl.ExpressionEvaluator.{ECMAScriptExpression, ECMAScriptFunction}
import shapeless.Poly1
import wom.values.WomValue

object EvaluateExpression extends Poly1 {
  implicit def script: Case.Aux[ECMAScriptExpression, (ParameterContext, Vector[ECMAScriptFunction]) => ErrorOr[WomValue]] = at[ECMAScriptExpression] { e =>
    (parameterContext: ParameterContext, expressionLib: Vector[ECMAScriptFunction]) =>
      ExpressionEvaluator.evalExpression(e, parameterContext, expressionLib)
  }

  implicit def function: Case.Aux[ECMAScriptFunction, (ParameterContext, Vector[ECMAScriptFunction]) => ErrorOr[WomValue]] = at[ECMAScriptFunction] { f =>
    (parameterContext: ParameterContext, expressionLib: Vector[ECMAScriptFunction]) =>
      ExpressionEvaluator.evalFunction(f, parameterContext, expressionLib)
  }
}
