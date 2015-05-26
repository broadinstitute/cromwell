package cromwell.binding.values

import cromwell.binding.WdlExpressionException
import cromwell.binding.types.WdlFloatType

import scala.util.{Try, Success, Failure}

case class WdlFloat(value: Double) extends WdlPrimitive {
  val wdlType = WdlFloatType
  override def add(rhs: WdlValue): Try[WdlValue] = {
    rhs match {
      case r:WdlFloat => Success(WdlFloat(value + r.value))
      case r:WdlInteger => Success(WdlFloat(value + r.value))
      case r:WdlString => Success(WdlString(value + r.value))
      case _ => invalid(s"$this + $rhs")
    }
  }
  override def subtract(rhs: WdlValue): Try[WdlValue] = {
    rhs match {
      case r:WdlFloat => Success(WdlFloat(value - r.value))
      case r:WdlInteger => Success(WdlFloat(value - r.value))
      case _ => invalid(s"$this - $rhs")
    }
  }
  override def multiply(rhs: WdlValue): Try[WdlValue] = {
    rhs match {
      case r:WdlFloat => Success(WdlFloat(value * r.value))
      case r:WdlInteger => Success(WdlFloat(value * r.value))
      case _ => invalid(s"$this * $rhs")
    }
  }
  override def divide(rhs: WdlValue): Try[WdlValue] = {
    rhs match {
      case r:WdlFloat if r.value == 0.0 => Failure(new WdlExpressionException("Divide by zero"))
      case r:WdlFloat => Success(WdlFloat(value / r.value))
      case r:WdlInteger if r.value == 0 => Failure(new WdlExpressionException("Divide by zero"))
      case r:WdlInteger => Success(WdlFloat(value / r.value))
      case _ => invalid(s"$this / $rhs")
    }
  }
  override def mod(rhs: WdlValue): Try[WdlValue] = {
    rhs match {
      case r:WdlFloat if r.value == 0.0 => Failure(new WdlExpressionException("Divide by zero"))
      case r:WdlFloat => Success(WdlFloat(value % r.value))
      case r:WdlInteger if r.value == 0 => Failure(new WdlExpressionException("Divide by zero"))
      case r:WdlInteger => Success(WdlFloat(value % r.value))
      case _ => invalid(s"$this % $rhs")
    }
  }
  override def equals(rhs: WdlValue): Try[WdlBoolean] = {
    rhs match {
      case r:WdlFloat => Success(WdlBoolean(value == r.value))
      case r:WdlInteger => Success(WdlBoolean(value == r.value))
      case _ => invalid(s"$this == $rhs")
    }
  }
  override def lessThan(rhs: WdlValue): Try[WdlBoolean] = {
    rhs match {
      case r:WdlFloat => Success(WdlBoolean(value < r.value))
      case r:WdlInteger => Success(WdlBoolean(value < r.value))
      case _ => invalid(s"$this < $rhs")
    }
  }
  override def greaterThan(rhs: WdlValue): Try[WdlBoolean] = {
    rhs match {
      case r:WdlFloat => Success(WdlBoolean(value > r.value))
      case r:WdlInteger => Success(WdlBoolean(value > r.value))
      case _ => invalid(s"$this > $rhs")
    }
  }
  override def unaryPlus: Try[WdlValue] = Success(WdlFloat(math.abs(value)))
  override def unaryMinus: Try[WdlValue] = Success(WdlFloat(-value))
  override def asString = value.toString
}
