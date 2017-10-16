package wom.values

import wom.types.WomBooleanType

import scala.util.{Success, Try}

object WomBoolean {
  val True = new WomBoolean(true)
  val False = new WomBoolean(false)
  def apply(value: Boolean) = if (value) WomBoolean.True else WomBoolean.False
  def unapply(b: WomBoolean): Option[Boolean] = Option(b.value)
}

/** The constructor is private to force access through the companion
  * object `apply` which ensures the use of one of the canonical instances.
  */
class WomBoolean private(val value: Boolean) extends WomPrimitive {
  val womType = WomBooleanType

  override def equals(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r:WomBoolean => Success(WomBoolean(value == r.value))
    case r: WomOptionalValue => evaluateIfDefined("==", r, equals)
    case _ => invalid(s"$value || $rhs")
  }

  override def lessThan(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r:WomBoolean => Success(WomBoolean(value < r.value))
    case r: WomOptionalValue => evaluateIfDefined("<", r, lessThan)
    case _ => invalid(s"$value < $rhs")
  }

  override def greaterThan(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r:WomBoolean => Success(WomBoolean(value > r.value))
    case r: WomOptionalValue => evaluateIfDefined(">", r, greaterThan)
    case _ => invalid(s"$value > $rhs")
  }

  override def or(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r:WomBoolean => Success(WomBoolean(value || r.value))
    case r: WomOptionalValue => evaluateIfDefined("||", r, or)
    case _ => invalid(s"$value || $rhs")
  }

  override def and(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r:WomBoolean => Success(WomBoolean(value && r.value))
    case r: WomOptionalValue => evaluateIfDefined("&&", r, and)
    case _ => invalid(s"$value && $rhs")
  }

  override def not: Try[WomValue] = Success(WomBoolean(!value))
  override def toWomString = value.toString
  override def toString = s"WomBoolean($value)"
}
