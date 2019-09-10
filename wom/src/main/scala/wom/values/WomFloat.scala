package wom.values

import wom.WomExpressionException
import wom.types.WomFloatType

import scala.util.{Failure, Success, Try}

case class WomFloat(value: BigDecimal) extends WomPrimitive {
  val womType = WomFloatType
  override def add(rhs: WomValue): Try[WomValue] = {
    rhs match {
      case r:WomFloat => Success(WomFloat(value + r.value))
      case r:WomInteger => Success(WomFloat(value + r.value))
      case r:WomString => Success(WomString(value + r.value))
      case r: WomOptionalValue => evaluateIfDefined("+", r, add)
      case _ => invalid(s"$this + $rhs")
    }
  }
  override def subtract(rhs: WomValue): Try[WomValue] = {
    rhs match {
      case r:WomFloat => Success(WomFloat(value - r.value))
      case r:WomInteger => Success(WomFloat(value - r.value))
      case r: WomOptionalValue => evaluateIfDefined("-", r, subtract)
      case _ => invalid(s"$this - $rhs")
    }
  }
  override def multiply(rhs: WomValue): Try[WomValue] = {
    rhs match {
      case r:WomFloat => Success(WomFloat(value * r.value))
      case r:WomInteger => Success(WomFloat(value * r.value))
      case r: WomOptionalValue => evaluateIfDefined("*", r, multiply)
      case _ => invalid(s"$this * $rhs")
    }
  }
  override def divide(rhs: WomValue): Try[WomValue] = {
    rhs match {
      case r:WomFloat if r.value == 0.0 => Failure(new WomExpressionException("Divide by zero"))
      case r:WomFloat => Success(WomFloat(value / r.value))
      case r:WomInteger if r.value == 0 => Failure(new WomExpressionException("Divide by zero"))
      case r:WomInteger => Success(WomFloat(value / r.value))
      case r: WomOptionalValue => evaluateIfDefined("/", r, divide)
      case _ => invalid(s"$this / $rhs")
    }
  }
  override def mod(rhs: WomValue): Try[WomValue] = {
    rhs match {
      case r:WomFloat if r.value == 0.0 => Failure(new WomExpressionException("Divide by zero"))
      case r:WomFloat => Success(WomFloat(value % r.value))
      case r:WomInteger if r.value == 0 => Failure(new WomExpressionException("Divide by zero"))
      case r:WomInteger => Success(WomFloat(value % r.value))
      case r: WomOptionalValue => evaluateIfDefined("%", r, mod)
      case _ => invalid(s"$this % $rhs")
    }
  }
  override def equals(rhs: WomValue): Try[WomBoolean] = {
    rhs match {
      case r:WomFloat => Success(WomBoolean(value == r.value))
      case r:WomInteger => Success(WomBoolean(value == r.value))
      case r: WomOptionalValue => evaluateIfDefined("==", r, equals)
      case _ => invalid(s"$this == $rhs")
    }
  }
  override def lessThan(rhs: WomValue): Try[WomBoolean] = {
    rhs match {
      case r:WomFloat => Success(WomBoolean(value < r.value))
      case r:WomInteger => Success(WomBoolean(value < r.value))
      case r: WomOptionalValue => evaluateIfDefined("<", r, lessThan)
      case _ => invalid(s"$this < $rhs")
    }
  }
  override def greaterThan(rhs: WomValue): Try[WomBoolean] = {
    rhs match {
      case r:WomFloat => Success(WomBoolean(value > r.value))
      case r:WomInteger => Success(WomBoolean(value > r.value))
      case r: WomOptionalValue => evaluateIfDefined(">", r, greaterThan)
      case _ => invalid(s"$this > $rhs")
    }
  }
  override def unaryPlus: Try[WomValue] = Success(WomFloat(value.abs))
  override def unaryMinus: Try[WomValue] = Success(WomFloat(-value))
  override def toWomString = value.toString
}
