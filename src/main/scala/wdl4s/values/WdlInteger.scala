package wdl4s.values

import wdl4s.WdlExpressionException
import wdl4s.types.WdlIntegerType

import scala.util.{Failure, Success, Try}

case class WdlInteger(value: Integer) extends WdlPrimitive {
  val wdlType = WdlIntegerType

  override def add(rhs: WdlValue): Try[WdlValue] = rhs match {
    case r:WdlInteger => Success(WdlInteger(value + r.value))
    case r:WdlString => Success(WdlString(value + r.value))
    case r:WdlFloat => Success(WdlFloat(value + r.value))
    case r: WdlOptionalValue => evaluateIfDefined(r, add)
    case _ => invalid(s"$value + $rhs")
  }

  override def subtract(rhs: WdlValue): Try[WdlValue] = rhs match {
    case r:WdlInteger => Success(WdlInteger(value - r.value))
    case r:WdlFloat => Success(WdlFloat(value - r.value))
    case r: WdlOptionalValue => evaluateIfDefined(r, subtract)
    case _ => invalid(s"$value - $rhs")
  }

  override def multiply(rhs: WdlValue): Try[WdlValue] = rhs match {
    case r:WdlInteger => Success(WdlInteger(value * r.value))
    case r:WdlFloat => Success(WdlFloat(value * r.value))
    case r: WdlOptionalValue => evaluateIfDefined(r, multiply)
    case _ => invalid(s"$value * $rhs")
  }

  override def divide(rhs: WdlValue): Try[WdlValue] = rhs match {
    case r:WdlInteger if r.value == 0 => Failure(new WdlExpressionException(s"Divide by zero error: $value / $rhs"))
    case r:WdlInteger => Success(WdlInteger(value / r.value))
    case r:WdlFloat if r.value == 0.toDouble => Failure(new WdlExpressionException(s"Divide by zero error: $value / $rhs"))
    case r:WdlFloat => Success(WdlFloat(value / r.value))
    case r: WdlOptionalValue => evaluateIfDefined(r, divide)
    case _ => invalid(s"$value / $rhs")
  }

  override def mod(rhs: WdlValue): Try[WdlValue] = rhs match {
    case r:WdlInteger if r.value == 0 => Failure(new WdlExpressionException(s"Divide by zero error: $value / $rhs"))
    case r:WdlInteger => Success(WdlInteger(value % r.value))
    case r:WdlFloat if r.value == 0.toDouble => Failure(new WdlExpressionException(s"Divide by zero error: $value / $rhs"))
    case r:WdlFloat => Success(WdlFloat(value % r.value))
    case r: WdlOptionalValue => evaluateIfDefined(r, mod)
    case _ => invalid(s"$value % $rhs")
  }

  override def equals(rhs: WdlValue): Try[WdlBoolean] = rhs match {
    case r:WdlInteger => Success(WdlBoolean(value == r.value))
    case r:WdlFloat => Success(WdlBoolean(value == r.value))
    case r: WdlOptionalValue => evaluateIfDefined(r, equals)
    case _ => invalid(s"$value == $rhs")
  }

  override def lessThan(rhs: WdlValue): Try[WdlBoolean] = rhs match {
    case r:WdlInteger => Success(WdlBoolean(value < r.value))
    case r:WdlFloat => Success(WdlBoolean(value < r.value))
    case r: WdlOptionalValue => evaluateIfDefined(r, lessThan)
    case _ => invalid(s"$value < $rhs")
  }

  override def greaterThan(rhs: WdlValue): Try[WdlBoolean] = rhs match {
    case r:WdlInteger => Success(WdlBoolean(value > r.value))
    case r:WdlFloat => Success(WdlBoolean(value > r.value))
    case r: WdlOptionalValue => evaluateIfDefined(r, greaterThan)
    case _ => invalid(s"$value > $rhs")
  }

  override def unaryPlus: Try[WdlValue] = Success(WdlInteger(math.abs(value)))
  override def unaryMinus: Try[WdlValue] = Success(WdlInteger(-value))
  override def toWdlString = value.toString
}
