package cwl

import common.validation.ErrorOr.ErrorOr
import cwl.ExpressionEvaluator.{ECMAScriptExpression, ECMAScriptFunction}
import shapeless.Poly1
import wom.values.WomValue

object EvaluateExpression extends Poly1 {
  implicit def script: Case.Aux[ECMAScriptExpression, (ParameterContext, ExpressionLib) => ErrorOr[WomValue]] = at[ECMAScriptExpression] { e =>
    (parameterContext: ParameterContext, expressionLib: ExpressionLib) =>
      ExpressionEvaluator.evalExpression(e, parameterContext, expressionLib)
  }

  implicit def function: Case.Aux[ECMAScriptFunction, (ParameterContext, ExpressionLib) => ErrorOr[WomValue]] = at[ECMAScriptFunction] { f =>
    (parameterContext: ParameterContext, expressionLib: ExpressionLib) =>
      ExpressionEvaluator.evalFunction(f, parameterContext, expressionLib)
  }
}
