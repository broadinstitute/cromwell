package cromwell.languages.util

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import shapeless.Poly1
import wom.expression.{IoFunctionSet, WomExpression}
import wom.values.WomValue

object ResolvedExecutableInputsPoly extends Poly1 {
  implicit def fromWdlValue: Case.Aux[WomValue, IoFunctionSet => ErrorOr[WomValue]] = at[WomValue] { wdlValue =>
    _: IoFunctionSet => wdlValue.validNel 
  }
  implicit def fromWomExpression: Case.Aux[WomExpression, IoFunctionSet => ErrorOr[WomValue]] = at[WomExpression] { womExpression =>
    ioFunctions: IoFunctionSet => womExpression.evaluateValue(Map.empty, ioFunctions).leftMap { errors => errors.map { e => s"Unable to evaluate expression '${womExpression.sourceString}': '$e'"}
    }
  }
}
