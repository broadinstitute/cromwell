package wdl.draft3.transforms.wdlom2wom.expression.values

import scala.language.implicitConversions
import common.validation.ErrorOr.ErrorOr
import simulacrum.typeclass
import wdl.draft3.transforms.wdlom2wom.expression.LinkedConsumedValue
import wdl.model.draft3.elements.ExpressionElement
import wom.values.WomValue

@typeclass
trait ValueEvaluator[A <: ExpressionElement] {
  def evaluateValue(a: A, inputs: Map[String, WomValue], linkedValues: Set[LinkedConsumedValue]): ErrorOr[WomValue]
}
