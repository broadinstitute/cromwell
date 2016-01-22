package wdl4s.values

import wdl4s.types.WdlBooleanType

import scala.util.{Success, Try}

object WdlBoolean {
  val True = new WdlBoolean(true)
  val False = new WdlBoolean(false)
  def apply(value: Boolean) = if (value) WdlBoolean.True else WdlBoolean.False
  def unapply(b: WdlBoolean): Option[Boolean] = Option(b.value)
}

/** The constructor is private to force access through the companion
  * object `apply` which ensures the use of one of the canonical instances.
  */
class WdlBoolean private(val value: Boolean) extends WdlPrimitive {
  val wdlType = WdlBooleanType

  override def equals(rhs: WdlValue): Try[WdlBoolean] = rhs match {
    case r:WdlBoolean => Success(WdlBoolean(value == r.value))
    case _ => invalid(s"$value || $rhs")
  }

  override def lessThan(rhs: WdlValue): Try[WdlBoolean] = rhs match {
    case r:WdlBoolean => Success(WdlBoolean(value < r.value))
    case _ => invalid(s"$value < $rhs")
  }

  override def greaterThan(rhs: WdlValue): Try[WdlBoolean] = rhs match {
    case r:WdlBoolean => Success(WdlBoolean(value > r.value))
    case _ => invalid(s"$value > $rhs")
  }

  override def or(rhs: WdlValue): Try[WdlBoolean] = rhs match {
    case r:WdlBoolean => Success(WdlBoolean(value || r.value))
    case _ => invalid(s"$value || $rhs")
  }

  override def and(rhs: WdlValue): Try[WdlBoolean] = rhs match {
    case r:WdlBoolean => Success(WdlBoolean(value && r.value))
    case _ => invalid(s"$value && $rhs")
  }

  override def not: Try[WdlValue] = Success(WdlBoolean(!value))
  override def toWdlString = value.toString
}
