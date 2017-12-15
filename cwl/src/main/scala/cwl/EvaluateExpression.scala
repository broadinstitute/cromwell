package cwl

import cwl.ExpressionEvaluator.{ECMAScriptExpression, ECMAScriptFunction}
import shapeless.Poly1
import wom.values.WomValue

import scala.util.Try

object EvaluateExpression extends Poly1 {
  implicit def script: Case.Aux[ECMAScriptExpression, ParameterContext => Try[WomValue]] = at[ECMAScriptExpression] { e =>
    (parameterContext: ParameterContext) =>
      ExpressionEvaluator.evalExpression(e, parameterContext)
  }

  implicit def function: Case.Aux[ECMAScriptFunction, ParameterContext => Try[WomValue]] = at[ECMAScriptFunction] { f =>
    (parameterContext: ParameterContext) =>
      ExpressionEvaluator.evalFunction(f, parameterContext)
  }
}
