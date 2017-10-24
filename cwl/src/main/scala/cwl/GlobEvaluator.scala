package cwl

import shapeless._
import wom.expression.IoFunctionSet
import wom.types.{WomArrayType, WomStringType}

import scala.Function._
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

  private type GlobHandler = ParameterContext => Seq[String]

  type Glob[A] = String

  def globPaths(glob: CommandOutputBinding.Glob,
                parameterContext: ParameterContext,
                ioFunctionSet: IoFunctionSet): Seq[String] = {
    val globs: Seq[String] = glob.fold(GlobToPaths).apply(parameterContext)
    getPaths(globs, ioFunctionSet)
  }

  /**
    * Find files relative to the output directory, using POSIX glob(3) pathname matching.
    *
    * @param globs The globs.
    * @return The paths that match the expression.
    */
  private def getPaths(globs: Seq[String], ioFunctionSet: IoFunctionSet): Seq[String] = {
    globs flatMap { glob =>
      //ioFunctionSet.glob(glob, "TODO: WOM: glob pattern here?")
      Vector(glob)
    }
  }

  object GlobToPaths extends Poly1 {
    implicit def caseECMAScript: Case.Aux[Expression, GlobHandler] = {
      at[Expression] { ecmaScript =>
        (parameterContext: ParameterContext) => {
          ecmaScript.fold(EvaluateExpression).apply(parameterContext) match {
            case WomArray(_, values) if values.isEmpty => Vector.empty
            case WomString(value) => Vector(value)
            case WomArray(WomArrayType(WomStringType), values) => values.map(_.valueString)
            case womValue =>
              throw new RuntimeException(
                s"Unexpected expression result: $womValue while evaluating expression '$ecmaScript' using inputs '${parameterContext.inputs}'")
          }
        }
      }
    }

    implicit def caseArrayString: Case.Aux[Array[String], GlobHandler] = {
      at[Array[String]] { const(_) }
    }

    implicit def caseString: Case.Aux[String, GlobHandler] = {
      at[String] { string => const(Vector(string)) }
    }
  }

}
