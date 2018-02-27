package cwl

import cats.syntax.validated._
import common.validation.ErrorOr._
import cwl.ExpressionEvaluator.eval
import cwl.internal.{CwlEcmaScriptDecoder, EcmaScriptUtil, EcmaScriptEncoder}
import eu.timepit.refined.api.Refined
import eu.timepit.refined.string.MatchesRegex
import shapeless.Witness
import wom.callable.RuntimeEnvironment
import wom.values.{WomFloat, WomInteger, WomString, WomValue}

// http://www.commonwl.org/v1.0/CommandLineTool.html#Expressions
object ExpressionEvaluator {
  // A code fragment wrapped in the $(...) syntax must be evaluated as a ECMAScript expression.
  val ECMAScriptExpressionWitness = Witness("""(?s)\s*\$\((.*)\)\s*""")
  val ECMAScriptExpressionRegex = ECMAScriptExpressionWitness.value.r
  type MatchesECMAScript = MatchesRegex[ECMAScriptExpressionWitness.T]
  type ECMAScriptExpression = String Refined MatchesECMAScript

  // A code fragment wrapped in the ${...} syntax must be evaluated as a ECMAScript function body for an anonymous,
  // zero-argument function.
  val ECMAScriptFunctionWitness = Witness("""(?s)\s*\$\{(.*)\}\s*""")
  val ECMAScriptFunctionRegex = ECMAScriptFunctionWitness.value.r
  type MatchesECMAFunction = MatchesRegex[ECMAScriptFunctionWitness.T]
  type ECMAScriptFunction = String Refined MatchesECMAFunction

  // This is not ECMAScript, just what CWL uses for interpolated expressions. The pattern is 'before$(expression)after'.
  // The regex uses a non-greedy quantifier on the 'before' to allow the expression to be processed from left to right,
  // and there are capturing groups around each portion. The official specification for interpolated strings is in the
  // last part of this section:
  // http://www.commonwl.org/v1.0/CommandLineTool.html#Parameter_references
  val InterpolatedStringWitness = Witness("""(?s)(.*?)\$\(([^\)]+)\)(.*)""")
  val InterpolatedStringRegex = InterpolatedStringWitness.value.r
  type MatchesInterpolatedString = MatchesRegex[InterpolatedStringWitness.T]
  type InterpolatedString = String Refined MatchesInterpolatedString

  def evalExpression(expression: ECMAScriptExpression, parameterContext: ParameterContext, expressionLib: ExpressionLib): ErrorOr[WomValue] = {
    expression.value match {
      case ECMAScriptExpressionRegex(script) =>
        // Nashorn doesn't like an expression floating around. So assign it to a variable and return that variable.
        val variableExpression =
          s"""|var expression_result = EXPRESSION_BODY;
              |expression_result
              |""".stripMargin.replace("EXPRESSION_BODY", script)
        
        eval(expressionFromParts(expressionLib, variableExpression), parameterContext)
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

  def evalInterpolatedString(string: InterpolatedString, parameterContext: ParameterContext, expressionLib: ExpressionLib): ErrorOr[WomValue] = {
    def interpolate(remaining: String, acc: String = ""): ErrorOr[String] = {
      remaining match {
        // The match is non-greedy in `before` so `expr` will always contain the first expression. e.g.:
        //
        // "foo $(bar) baz $(qux) quux" would match with before = "foo ", expr = "bar", after = " baz $(qux) quux"
        case InterpolatedStringRegex(before, expr, after) =>
          eval(expressionFromParts(expressionLib, expr), parameterContext) flatMap { v =>
            interpolate(after, acc + before + v.valueString)
          }
        case r => (acc + r).validNel
      }
    }

    string.value match {
      case InterpolatedStringRegex(_, _, _) => interpolate(string.value) map WomString.apply

      case unmatched =>
        s"Expression '$unmatched' was unable to be matched to regex '${InterpolatedStringWitness.value}'".invalidNel
    }
  }

  private lazy val cwlJsEncoder = new EcmaScriptEncoder()
  private lazy val cwlJsDecoder = new CwlEcmaScriptDecoder()

  def eval(expr: String, parameterContext: ParameterContext): ErrorOr[WomValue] = {
    val (rawValues, mapValues) = paramValues(parameterContext)
    EcmaScriptUtil.evalStructish(expr, rawValues, mapValues, cwlJsEncoder, cwlJsDecoder)
  }

  def eval(expr: Expression, parameterContext: ParameterContext, expressionLib: ExpressionLib): ErrorOr[WomValue] = {
    expr.fold(EvaluateExpression).apply(parameterContext, expressionLib)
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
      "cores" -> WomInteger(runtime.cores),
      "ram" -> WomFloat(runtime.ram),
      "outdirSize" -> WomFloat(runtime.outputPathSize.toDouble),
      "tmpdirSize" -> WomFloat(runtime.tempPathSize.toDouble)
    )
  }
}
