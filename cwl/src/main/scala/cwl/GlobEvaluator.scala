package cwl

import cats.syntax.validated._
import common.validation.ErrorOr._
import shapeless._
import wom.types.{WomArrayType, WomStringType}
import wom.values._

/*
CommandOutputBinding.glob:
Find files relative to the output directory, using POSIX glob(3) pathname matching. If an array is provided, find
files that match any pattern in the array. If an expression is provided, the expression must return a string or an
array of strings, which will then be evaluated as one or more glob patterns. Must only match and return files which
actually exist.

http://www.commonwl.org/v1.0/CommandLineTool.html#CommandOutputBinding
 */
object GlobEvaluator {

  type GlobEvaluatorFunction = (ParameterContext, ExpressionLib) => ErrorOr[List[String]]

  def globs(globOption: Option[Glob], parameterContext: ParameterContext, expressionLib: ExpressionLib): ErrorOr[List[String]] = {
    globOption map {
      globs(_, parameterContext, expressionLib)
    } getOrElse {
      Nil.valid
    }
  }

  def globs(glob: Glob, parameterContext: ParameterContext, expressionLib: ExpressionLib): ErrorOr[List[String]] = {
    glob.fold(GlobEvaluator.GlobEvaluatorPoly).apply(parameterContext, expressionLib)
  }

  object GlobEvaluatorPoly extends Poly1 {
    implicit def caseStringOrExpression: Case.Aux[StringOrExpression, GlobEvaluatorFunction] = {
      at {
        _.fold(this)
      }
    }

    implicit def caseExpression: Case.Aux[Expression, GlobEvaluatorFunction] = {
      at {
        expression =>
          (parameterContext, expressionLib) => {
            ExpressionEvaluator.eval(expression, parameterContext) flatMap {
              case WomArray(_, values) if values.isEmpty => Nil.valid
              case WomString(value) => List(value).valid
              case WomArray(WomArrayType(WomStringType), values) => values.toList.map(_.valueString).valid
              case womValue =>
                val message = s"Unexpected expression result: $womValue while evaluating expression '$expression' " +
                  s"using inputs '$parameterContext'"
                message.invalidNel
            }
          }
      }
    }

    implicit def caseArrayString: Case.Aux[Array[String], GlobEvaluatorFunction] = {
      at { array => (_, _) => array.toList.valid }
    }

    implicit def caseString: Case.Aux[String, GlobEvaluatorFunction] = {
      at { string => (_, _) => List(string).valid }
    }
  }

}
