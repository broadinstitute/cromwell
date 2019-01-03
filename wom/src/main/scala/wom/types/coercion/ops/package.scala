package wom.types.coercion

import common.validation.ErrorOr.ErrorOr

package object ops {
  implicit class AnyCoercer(val a: Any) extends AnyVal {
    def coercionDefined[A](implicit coercer: WomTypeCoercer[A]): Boolean = coercer.coercionDefined(a)
    def coerceToType[A](implicit coercer: WomTypeCoercer[A]): ErrorOr[A] = coercer.coerceToType(a)
  }
}
