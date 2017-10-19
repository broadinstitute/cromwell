package wom.values

import wom.WdlExpressionException
import wom.types.WdlFloatType

import scala.util.{Failure, Success, Try}

case class WdlFloat(value: Double) extends WdlPrimitive {
  val wdlType = WdlFloatType
  override def add(rhs: WdlValue): Try[WdlValue] = {
    rhs match {
      case r:WdlFloat => Success(WdlFloat(value + r.value))
      case r:WdlInteger => Success(WdlFloat(value + r.value))
      case r:WdlString => Success(WdlString(value + r.value))
      case r: WdlOptionalValue => evaluateIfDefined("+", r, add)
      case _ => invalid(s"$this + $rhs")
    }
  }
  override def subtract(rhs: WdlValue): Try[WdlValue] = {
    rhs match {
      case r:WdlFloat => Success(WdlFloat(value - r.value))
      case r:WdlInteger => Success(WdlFloat(value - r.value))
      case r: WdlOptionalValue => evaluateIfDefined("-", r, subtract)
      case _ => invalid(s"$this - $rhs")
    }
  }
  override def multiply(rhs: WdlValue): Try[WdlValue] = {
    rhs match {
      case r:WdlFloat => Success(WdlFloat(value * r.value))
      case r:WdlInteger => Success(WdlFloat(value * r.value))
      case r: WdlOptionalValue => evaluateIfDefined("*", r, multiply)
      case _ => invalid(s"$this * $rhs")
    }
  }
  override def divide(rhs: WdlValue): Try[WdlValue] = {
    rhs match {
      case r:WdlFloat if r.value == 0.0 => Failure(new WdlExpressionException("Divide by zero"))
      case r:WdlFloat => Success(WdlFloat(value / r.value))
      case r:WdlInteger if r.value == 0 => Failure(new WdlExpressionException("Divide by zero"))
      case r:WdlInteger => Success(WdlFloat(value / r.value))
      case r: WdlOptionalValue => evaluateIfDefined("/", r, divide)
      case _ => invalid(s"$this / $rhs")
    }
  }
  override def mod(rhs: WdlValue): Try[WdlValue] = {
    rhs match {
      case r:WdlFloat if r.value == 0.0 => Failure(new WdlExpressionException("Divide by zero"))
      case r:WdlFloat => Success(WdlFloat(value % r.value))
      case r:WdlInteger if r.value == 0 => Failure(new WdlExpressionException("Divide by zero"))
      case r:WdlInteger => Success(WdlFloat(value % r.value))
      case r: WdlOptionalValue => evaluateIfDefined("%", r, mod)
      case _ => invalid(s"$this % $rhs")
    }
  }
  override def equals(rhs: WdlValue): Try[WdlBoolean] = {
    rhs match {
      case r:WdlFloat => Success(WdlBoolean(value == r.value))
      case r:WdlInteger => Success(WdlBoolean(value == r.value))
      case r: WdlOptionalValue => evaluateIfDefined("==", r, equals)
      case _ => invalid(s"$this == $rhs")
    }
  }
  override def lessThan(rhs: WdlValue): Try[WdlBoolean] = {
    rhs match {
      case r:WdlFloat => Success(WdlBoolean(value < r.value))
      case r:WdlInteger => Success(WdlBoolean(value < r.value))
      case r: WdlOptionalValue => evaluateIfDefined("<", r, lessThan)
      case _ => invalid(s"$this < $rhs")
    }
  }
  override def greaterThan(rhs: WdlValue): Try[WdlBoolean] = {
    rhs match {
      case r:WdlFloat => Success(WdlBoolean(value > r.value))
      case r:WdlInteger => Success(WdlBoolean(value > r.value))
      case r: WdlOptionalValue => evaluateIfDefined(">", r, greaterThan)
      case _ => invalid(s"$this > $rhs")
    }
  }
  override def unaryPlus: Try[WdlValue] = Success(WdlFloat(math.abs(value)))
  override def unaryMinus: Try[WdlValue] = Success(WdlFloat(-value))
  override def toWdlString = value.toString
}
