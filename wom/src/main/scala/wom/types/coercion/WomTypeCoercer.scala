package wom.types.coercion

import common.validation.ErrorOr.ErrorOr

trait WomTypeCoercer[A] {
  def toDisplayString: String
  def coercionDefined(any: Any): Boolean
  def coerceToType(any: Any): ErrorOr[A]
}
