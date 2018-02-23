package wdl.draft3.transforms.wdlom2wom.expression.types

import scala.language.implicitConversions

import common.validation.ErrorOr.ErrorOr
import simulacrum.typeclass
import wdl.draft3.transforms.wdlom2wom.expression.LinkedConsumedValue
import wdl.model.draft3.elements.ExpressionElement
import wom.types.WomType

@typeclass
trait TypeEvaluator[A <: ExpressionElement] {
  def evaluateType(a: A, linkedValues: Set[LinkedConsumedValue]): ErrorOr[WomType]
}
