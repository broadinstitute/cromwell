package wom.values

import wom.WomExpressionException
import wom.types.WomBigDecimalType

import scala.util.{Failure, Success, Try}

case class WomBigDecimal(value: BigDecimal) extends WomPrimitive {
  val womType = WomBigDecimalType

  override def add(rhs: WomValue): Try[WomValue] = rhs match {
    case r: WomInteger => Success(WomBigDecimal(value + r.value))
    case r: WomBigDecimal => Success(WomBigDecimal(value + r.value))
    case r: WomString => Success(WomString(value + r.value))
    case r: WomFloat => Success(WomBigDecimal(value + r.value))
    case r: WomOptionalValue => evaluateIfDefined("+", r, add)
    case _ => invalid(s"$value + $rhs")
  }

  override def subtract(rhs: WomValue): Try[WomValue] = rhs match {
    case r: WomInteger => Success(WomBigDecimal(value - r.value))
    case r: WomBigDecimal => Success(WomBigDecimal(value - r.value))
    case r: WomFloat => Success(WomBigDecimal(value - r.value))
    case r: WomOptionalValue => evaluateIfDefined("-", r, subtract)
    case _ => invalid(s"$value - $rhs")
  }

  override def multiply(rhs: WomValue): Try[WomValue] = rhs match {
    case r: WomInteger => Success(WomBigDecimal(value * r.value))
    case r: WomBigDecimal => Success(WomBigDecimal(value * r.value))
    case r: WomFloat => Success(WomBigDecimal(value * r.value))
    case r: WomOptionalValue => evaluateIfDefined("*", r, multiply)
    case _ => invalid(s"$value * $rhs")
  }

  override def divide(rhs: WomValue): Try[WomValue] = rhs match {
    case r: WomInteger if r.value == 0 => Failure(new WomExpressionException(s"Divide by zero error: $value / $rhs"))
    case r: WomInteger => Success(WomBigDecimal(value / r.value))
    case r: WomBigDecimal if r.value == 0 => Failure(new WomExpressionException(s"Divide by zero error: $value / $rhs"))
    case r: WomBigDecimal => Success(WomBigDecimal(value / r.value))
    case r: WomFloat if r.value == 0.toDouble => Failure(new WomExpressionException(s"Divide by zero error: $value / $rhs"))
    case r: WomFloat => Success(WomBigDecimal(value / r.value))
    case r: WomOptionalValue => evaluateIfDefined("/", r, divide)
    case _ => invalid(s"$value / $rhs")
  }

  override def mod(rhs: WomValue): Try[WomValue] = rhs match {
    case r: WomInteger if r.value == 0 => Failure(new WomExpressionException(s"Divide by zero error: $value / $rhs"))
    case r: WomInteger => Success(WomBigDecimal(value % r.value))
    case r: WomBigDecimal if r.value == 0 => Failure(new WomExpressionException(s"Divide by zero error: $value / $rhs"))
    case r: WomBigDecimal => Success(WomBigDecimal(value % r.value))
    case r: WomFloat if r.value == 0.toDouble => Failure(new WomExpressionException(s"Divide by zero error: $value / $rhs"))
    case r: WomFloat => Success(WomBigDecimal(value % r.value))
    case r: WomOptionalValue => evaluateIfDefined("%", r, mod)
    case _ => invalid(s"$value % $rhs")
  }

  override def equals(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r: WomInteger => Success(WomBoolean(value == r.value))
    case r: WomBigDecimal => Success(WomBoolean(value == r.value))
    case r: WomFloat => Success(WomBoolean(value == r.value))
    case r: WomOptionalValue => evaluateIfDefined("==", r, equals)
    case _ => invalid(s"$value == $rhs")
  }

  override def lessThan(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r: WomInteger => Success(WomBoolean(value < r.value))
    case r: WomBigDecimal => Success(WomBoolean(value < r.value))
    case r: WomFloat => Success(WomBoolean(value < r.value))
    case r: WomOptionalValue => evaluateIfDefined("<", r, lessThan)
    case _ => invalid(s"$value < $rhs")
  }

  override def greaterThan(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r: WomInteger => Success(WomBoolean(value > r.value))
    case r: WomBigDecimal => Success(WomBoolean(value > r.value))
    case r: WomFloat => Success(WomBoolean(value > r.value))
    case r: WomOptionalValue => evaluateIfDefined(">", r, greaterThan)
    case _ => invalid(s"$value > $rhs")
  }

  override def unaryPlus: Try[WomValue] = Success(WomBigDecimal(value.abs))

  override def unaryMinus: Try[WomValue] = Success(WomBigDecimal(-value))

  override def toWomString = value.toString
}
