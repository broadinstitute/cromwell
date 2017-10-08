package wdl4s.cwl

import cwl.ExpressionEvaluator.{ECMAScriptExpression, ECMAScriptFunction}
import cwl.{ExpressionEvaluator, ParameterContext}
import shapeless.Poly1

object EvaluateExpression extends Poly1 {
  implicit def script = at[ECMAScriptExpression] { e =>
    (parameterContext: ParameterContext) =>
      ExpressionEvaluator.evalExpression(e, parameterContext)
  }

  implicit def function = at[ECMAScriptFunction] { f =>
    (parameterContext: ParameterContext) =>
      ExpressionEvaluator.evalFunction(f, parameterContext)
  }
}
