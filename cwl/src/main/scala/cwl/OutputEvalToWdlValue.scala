package cwl

import shapeless._
import wom.values.{WomString, WomValue}

import scala.Function._

object OutputEvalToWdlValue extends Poly1 {

  type OutputEvalHandler = ParameterContext => WomValue

  implicit def caseECMAScript: Case.Aux[Expression, OutputEvalHandler] = {
    at[Expression] { ecmaScript =>
      (parameterContext: ParameterContext) =>
        ecmaScript.fold(EvaluateExpression).apply(parameterContext)
    }
  }

  implicit def caseString: Case.Aux[String, OutputEvalHandler] = {
    at[String] { string => const(WomString(string)) }
  }
}


