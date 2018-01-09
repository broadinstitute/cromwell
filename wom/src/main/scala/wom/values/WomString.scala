package wom.values

import org.apache.commons.text.StringEscapeUtils
import wom.types.WomStringType

import scala.util.{Success, Try}

case class WomString(value: String) extends WomPrimitive {
  val womType = WomStringType

  override def add(rhs: WomValue): Try[WomValue] = rhs match {
    case r: WomString => Success(WomString(value + r.value))
    case r: WomInteger => Success(WomString(value + r.value))
    case r: WomFloat => Success(WomString(value + r.value))
    case r: WomPrimitiveFile => Success(WomString(value + r.value))
    case r: WomOptionalValue => evaluateIfDefined("+", r, add)
    case _ => invalid(s"$value + $rhs")
  }

  override def equals(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r: WomPrimitiveFile => Success(WomBoolean(value == r.value))
    case r: WomString => Success(WomBoolean(value == r.value))
    case r: WomOptionalValue => evaluateIfDefined("==", r, equals)
    case _ => invalid(s"$value == $rhs")
  }

  override def lessThan(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r: WomString => Success(WomBoolean(value < r.value))
    case r: WomOptionalValue => evaluateIfDefined("<", r, lessThan)
    case _ => invalid(s"$value < $rhs")
  }

  override def greaterThan(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r: WomString => Success(WomBoolean(value > r.value))
    case r: WomOptionalValue => evaluateIfDefined(">", r, greaterThan)
    case _ => invalid(s"$value > $rhs")
  }

  override def toWomString = "\"" + StringEscapeUtils.escapeJava(value) + "\""
  override def valueString = value
}
