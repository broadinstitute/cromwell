package cwl

import cats.syntax.validated._
import common.validation.ErrorOr._
import cwl.internal.{CwlEcmaScriptDecoder, EcmaScriptEncoder, EcmaScriptUtil}
import eu.timepit.refined.api.Refined
import eu.timepit.refined.string.MatchesRegex
import shapeless.Witness
import wom.callable.RuntimeEnvironment
import wom.values.{WomFloat, WomInteger, WomString, WomValue}

// http://www.commonwl.org/v1.0/CommandLineTool.html#Expressions
object ExpressionEvaluator {
  // A code fragment wrapped in the $(...) syntax must be evaluated as a ECMAScript expression.
  /*
  There are several places in cwltool where a simple contains check is used:
  - https://github.com/common-workflow-language/cwltool/blob/353dbed/cwltool/builder.py#L157
  - https://github.com/common-workflow-language/cwltool/blob/353dbed/cwltool/command_line_tool.py#L669
  - https://github.com/common-workflow-language/cwltool/blob/353dbed/cwltool/expression.py#L219
  - https://github.com/common-workflow-language/cwltool/blob/353dbed/cwltool/process.py#L645
   */
  val ECMAScriptExpressionWitness = Witness("""(?s).*\$\(.*""")
  val ECMAScriptExpressionRegex = ECMAScriptExpressionWitness.value.r
  type MatchesECMAScriptExpression = MatchesRegex[ECMAScriptExpressionWitness.T]
  type ECMAScriptExpression = String Refined MatchesECMAScriptExpression

  // A code fragment wrapped in the ${...} syntax must be evaluated as a ECMAScript function body for an anonymous,
  // zero-argument function.
  val ECMAScriptFunctionWitness = Witness("""(?s)\s*\$\{(.*)\}\s*""")
  val ECMAScriptFunctionRegex = ECMAScriptFunctionWitness.value.r
  type MatchesECMAScriptFunction = MatchesRegex[ECMAScriptFunctionWitness.T]
  type ECMAScriptFunction = String Refined MatchesECMAScriptFunction

  def eval(expr: Expression, parameterContext: ParameterContext): ErrorOr[WomValue] = {
    expr.fold(EvaluateExpression).apply(parameterContext)
  }

  def evalExpression(expression: ECMAScriptExpression)(parameterContext: ParameterContext): ErrorOr[WomValue] = {
    def evaluator(string: String): ErrorOr[WomValue] = {
      eval(string, parameterContext)
    }

    ExpressionInterpolator.interpolate(expression.value, evaluator, strip_whitespace = false)
  }

  def evalFunction(function: ECMAScriptFunction)(parameterContext: ParameterContext): ErrorOr[WomValue] = {
    function.value match {
      case ECMAScriptFunctionRegex(script) =>

        val functionExpression =
          s"""|(function() {
              |FUNCTION_BODY
              |})();
              |""".stripMargin.replace("FUNCTION_BODY", script)

        eval(functionExpression, parameterContext)
      case unmatched =>
        s"Expression '$unmatched' was unable to be matched to regex '${ECMAScriptFunctionWitness.value}'".invalidNel
    }
  }

  private lazy val cwlJsDecoder = new CwlEcmaScriptDecoder()

  def eval(expr: String, parameterContext: ParameterContext): ErrorOr[WomValue] = {
    val script = if (parameterContext.expressionLib.isEmpty) {
      expr
    } else {
      parameterContext.expressionLib.mkString("", ";", s";$expr")
    }
    val (rawVals, mapVals) = paramValues(parameterContext)
    EcmaScriptUtil.evalStructish(
      script,
      rawVals,
      mapVals,
      new EcmaScriptEncoder(),
      cwlJsDecoder
    )
  }

  def paramValues(parameterContext: ParameterContext): ((String, WomValue), Map[String, Map[String, WomValue]]) = {
    (
        "self" -> parameterContext.self
      ,
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
      "cores" -> WomInteger(runtime.cores.value),
      "ram" -> WomFloat(runtime.ram),
      "outdirSize" -> WomFloat(runtime.outputPathSize.toDouble),
      "tmpdirSize" -> WomFloat(runtime.tempPathSize.toDouble)
    )
  }
}
