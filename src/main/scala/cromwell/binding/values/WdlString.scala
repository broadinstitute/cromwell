package cromwell.binding.values

import cromwell.binding.types.WdlStringType

import scala.util.{Success, Try}

case class WdlString(value: String) extends WdlPrimitive {
  val wdlType = WdlStringType
  override def add(rhs: WdlValue): Try[WdlValue] = {
    rhs match {
      case r:WdlString => Success(WdlString(value + r.value))
      case r:WdlInteger => Success(WdlString(value + r.value))
      case r:WdlFloat => Success(WdlString(value + r.value))
      case _ => invalid(s"$value + $rhs")
    }
  }
  override def equals(rhs: WdlValue): Try[WdlBoolean] = {
    rhs match {
      case r:WdlString => Success(WdlBoolean(value == r.value))
      case _ => invalid(s"$value == $rhs")
    }
  }
  override def lessThan(rhs: WdlValue): Try[WdlBoolean] = {
    rhs match {
      case r:WdlString => Success(WdlBoolean(value < r.value))
      case _ => invalid(s"$value < $rhs")
    }
  }
  override def greaterThan(rhs: WdlValue): Try[WdlBoolean] = {
    rhs match {
      case r:WdlString => Success(WdlBoolean(value > r.value))
      case _ => invalid(s"$value > $rhs")
    }
  }
  override def asString = value
}
