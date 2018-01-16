package cwl

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import eu.timepit.refined.api.Refined
import eu.timepit.refined.string.MatchesRegex
import shapeless.Witness
import wom.callable.RuntimeEnvironment
import wom.util.JsUtil
import wom.values.{WomFloat, WomInteger, WomString, WomValue}

// http://www.commonwl.org/v1.0/CommandLineTool.html#Expressions
object ExpressionEvaluator {
  // A code fragment wrapped in the $(...) syntax must be evaluated as a ECMAScript expression.
  val ECMAScriptExpressionWitness = Witness("""(?s)\s*\$\((.*)\)\s*""")
  val ECMAScriptExpressionRegex = ECMAScriptExpressionWitness.value.r
  type MatchesECMAScript = MatchesRegex[ECMAScriptExpressionWitness.T]
  type ECMAScriptExpression = String Refined MatchesRegex[ECMAScriptExpressionWitness.T]

  // A code fragment wrapped in the ${...} syntax must be evaluated as a ECMAScript function body for an anonymous,
  // zero-argument function.
  val ECMAScriptFunctionWitness = Witness("""(?s)\s*\$\{(.*)\}\s*""")
  val ECMAScriptFunctionRegex = ECMAScriptFunctionWitness.value.r
  type ECMAScriptFunction = String Refined MatchesRegex[ECMAScriptFunctionWitness.T]
  type MatchesECMAFunction = MatchesRegex[ECMAScriptFunctionWitness.T]

  def evalExpression(expression: ECMAScriptExpression, parameterContext: ParameterContext, expressionLib: ExpressionLib): ErrorOr[WomValue] = {
    expression.value match {
      case ECMAScriptExpressionRegex(script) =>
        eval(expressionFromParts(expressionLib, script), parameterContext)
      case unmatched =>
        s"Expression '$unmatched' was unable to be matched to regex '${ECMAScriptExpressionWitness.value}'".invalidNel
    }
  }

  def expressionFromParts(lib: ExpressionLib, script: String) = {
    val expressionScript: String = lib.mkString(";")
    s"$expressionScript;$script"
  }

  def evalFunction(function: ECMAScriptFunction, parameterContext: ParameterContext, expressionLib: ExpressionLib): ErrorOr[WomValue] = {
    function.value match {
      case ECMAScriptFunctionRegex(script) =>

        val functionExpression =
          s"""|(function() {
              |FUNCTION_BODY
              |})();
              |""".stripMargin.replace("FUNCTION_BODY", script)

        eval(expressionFromParts(expressionLib, functionExpression), parameterContext)
      case unmatched =>
        s"Expression '$unmatched' was unable to be matched to regex '${ECMAScriptFunctionWitness.value}'".invalidNel
    }
  }

  private lazy val cwlJsEncoder = new CwlJsEncoder()
  private lazy val cwlJsDecoder = new CwlJsDecoder()

  def eval(expr: String, parameterContext: ParameterContext): ErrorOr[WomValue] = {
    val (rawValues, mapValues) = paramValues(parameterContext)
    JsUtil.evalStructish(expr, rawValues, mapValues, cwlJsEncoder, cwlJsDecoder)
  }

  def eval(expr: Expression, parameterContext: ParameterContext, expressionLib: ExpressionLib): ErrorOr[WomValue] = {
    expr.fold(EvaluateExpression).apply(parameterContext, expressionLib)
  }

  def paramValues(parameterContext: ParameterContext): (Map[String, WomValue], Map[String, Map[String, WomValue]]) = {
    (
      Map(
        "self" -> parameterContext.self
      ),
      Map(
        "inputs" -> parameterContext.inputs,
        "runtime" -> parameterContext.runtimeOption.map(cwlMap).getOrElse(Map.empty)
      )
    )
  }

  def cwlMap(runtime: RuntimeEnvironment): Map[String, WomValue] = {
    Map(
      "outdir" -> WomString(runtime.outputPath),
      "tmpdir" -> WomString(runtime.tempPath),
      "cores" -> WomInteger(runtime.cores),
      "ram" -> WomFloat(runtime.ram),
      "outdirSize" -> WomFloat(runtime.outputPathSize.toDouble),
      "tmpdirSize" -> WomFloat(runtime.tempPathSize.toDouble)
    )
  }
}
