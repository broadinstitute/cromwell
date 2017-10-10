package cromwell.engine.workflow.lifecycle

import cats.syntax.validated._
import lenthall.validation.ErrorOr.ErrorOr
import shapeless.Poly1
import wdl.values.WdlValue
import wom.expression.{IoFunctionSet, WomExpression}

object ResolvedExecutableInputsPoly extends Poly1 {
  implicit def fromWdlValue: Case.Aux[WdlValue, IoFunctionSet => ErrorOr[WdlValue]] = at[WdlValue] { wdlValue =>
    _: IoFunctionSet => wdlValue.validNel 
  }
  implicit def fromWomExpression: Case.Aux[WomExpression, IoFunctionSet => ErrorOr[WdlValue]] = at[WomExpression] { womExpression =>
    ioFunctions: IoFunctionSet => womExpression.evaluateValue(Map.empty, ioFunctions)
  }
}
