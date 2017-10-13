package cwl

import eu.timepit.refined.api.Refined
import eu.timepit.refined.string.MatchesRegex
import shapeless.Witness
import wom.util.JsUtil
import wom.values.WdlValue

// http://www.commonwl.org/v1.0/CommandLineTool.html#Expressions
object ExpressionEvaluator {
  // A code fragment wrapped in the $(...) syntax must be evaluated as a ECMAScript expression.
  val ECMAScriptExpressionWitness = Witness("""\$\((.*)\)""")
  type MatchesECMAScript = MatchesRegex[ECMAScriptExpressionWitness.T]
  type ECMAScriptExpression = String Refined MatchesRegex[ECMAScriptExpressionWitness.T]

  // A code fragment wrapped in the ${...} syntax must be evaluated as a ECMAScript function body for an anonymous,
  // zero-argument function.
  val ECMAScriptFunctionWitness = Witness("""\$\{(.*)\}""")
  type ECMAScriptFunction = String Refined MatchesRegex[ECMAScriptFunctionWitness.T]
  type MatchesECMAFunction = MatchesRegex[ECMAScriptFunctionWitness.T]


  def evalExpression(expression: ECMAScriptExpression, parameterContext: ParameterContext): WdlValue = {
    val ECMAScriptExpressionRegex = ECMAScriptExpressionWitness.value.r
    expression.value match {
      case ECMAScriptExpressionRegex(script) => JsUtil.eval(script, paramValues(parameterContext))
      case _ => throw new RuntimeException(s"Expression was unable to be matched to Regex. This is never supposed to happen thanks to our JSON parsing library")
    }
  }

  def evalFunction(function: ECMAScriptFunction, parameterContext: ParameterContext): WdlValue = {
    val ECMAScriptFunctionRegex = ECMAScriptFunctionWitness.value.r
    function.value match {
      case ECMAScriptFunctionRegex(script) =>
        val functionExpression =
          s"""|(function() {
              |FUNCTION_BODY
              |})();
              |""".stripMargin.replaceAll("FUNCTION_BODY", script)

        JsUtil.eval(functionExpression, paramValues(parameterContext))
      case _ => throw new RuntimeException(s"Expression was unable to be matched to Regex. This is never supposed to happen thanks to our JSON parsing library")
    }
  }

  def paramValues(parameterContext: ParameterContext) =
    Map(
      "inputs" -> parameterContext.inputs,
      "runtime" -> parameterContext.runtime,
      "self" -> parameterContext.self
    )

}
