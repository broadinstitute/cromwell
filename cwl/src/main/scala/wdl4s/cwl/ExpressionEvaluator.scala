package wdl4s.cwl

import wdl4s.wdl.util.JsUtil
import wdl4s.wdl.values.WdlValue

// http://www.commonwl.org/v1.0/CommandLineTool.html#Expressions
object ExpressionEvaluator {
  // A code fragment wrapped in the $(...) syntax must be evaluated as a ECMAScript expression.
  private val EcmaScriptRegex = """\$\((.*)\)""".r

  // A code fragment wrapped in the ${...} syntax must be evaluated as a ECMAScript function body for an anonymous,
  // zero-argument function.
  private val EcmaFunctionRegex = """\$\{(.*)\}""".r

  def evalExpression(expression: String, parameterContext: ParameterContext): WdlValue = {
    val evalValues = Map(
      "inputs" -> parameterContext.inputs,
      "runtime" -> parameterContext.runtime,
      "self" -> parameterContext.self
    )

    // TODO: WOM: Are we supposed to be trimming the expression just in case before matching?
    expression.trim match {
      case EcmaScriptRegex(ecmaScriptExpression) => JsUtil.eval(ecmaScriptExpression, evalValues)
      case EcmaFunctionRegex(functionBody) =>
        val functionExpression =
          s"""|(function() {
              |FUNCTION_BODY
              |})();
              |""".stripMargin.replaceAll("FUNCTION_BODY", functionBody)
        JsUtil.eval(functionExpression, evalValues)
      case _ => throw new RuntimeException(s"TODO: WOM: Got unexpected expression '$expression'")
    }
  }
}
