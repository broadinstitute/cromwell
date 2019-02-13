package wom.types

import spray.json.{JsNumber, JsString}
import wom.values.{WomFloat, WomInteger, WomLong, WomString}

import scala.util.{Success, Try}

case object WomFloatType extends WomPrimitiveType {
  val stableName: String = "Float"

  override protected def coercion = {
    case f: Float => WomFloat(f.toDouble)
    case d: Double => WomFloat(d)
    case i: Integer => WomFloat(i.toDouble)
    case n: JsNumber => WomFloat(n.value.doubleValue())
    case f: WomFloat => f
    case i: WomInteger => WomFloat(i.value.toDouble)
    case l: WomLong => WomFloat(l.value.toDouble)
    case s: String => WomFloat(s.toDouble)
    case s: JsString => WomFloat(s.value.toDouble)
    case s: WomString => WomFloat(s.value.toDouble)
  }

  private def binaryOperator(rhs: WomType, symbol: String): Try[WomType] = rhs match {
    case WomIntegerType => Success(WomFloatType)
    case WomFloatType => Success(WomFloatType)
    case WomOptionalType(memberType) => binaryOperator(memberType, symbol)
    case _ => invalid(s"$this $symbol $rhs")
  }

  private def comparisonOperator(rhs: WomType, symbol: String): Try[WomType] = rhs match {
    case wct:WomCoproductType => wct.typeExists(WomStringType)
    case WomIntegerType => Success(WomBooleanType)
    case WomFloatType => Success(WomBooleanType)
    case WomOptionalType(memberType) => comparisonOperator(memberType, symbol)
    case _ => invalid(s"$this $symbol $rhs")
  }

  override def add(rhs: WomType): Try[WomType] = rhs match {
    case WomStringType => Success(WomStringType)
    case WomOptionalType(memberType) => add(memberType)
    case t => binaryOperator(t, "+")
  }

  override def subtract(rhs: WomType): Try[WomType] = binaryOperator(rhs, "-")
  override def multiply(rhs: WomType): Try[WomType] = binaryOperator(rhs, "*")
  override def divide(rhs: WomType): Try[WomType] = binaryOperator(rhs, "/")
  override def mod(rhs: WomType): Try[WomType] = binaryOperator(rhs, "%")
  override def equalsType(rhs: WomType): Try[WomType] = comparisonOperator(rhs, "==")
  override def lessThan(rhs: WomType): Try[WomType] = comparisonOperator(rhs, "<")
  override def greaterThan(rhs: WomType): Try[WomType] = comparisonOperator(rhs, ">")
  override def unaryPlus: Try[WomType] = Success(WomFloatType)
  override def unaryMinus: Try[WomType] = Success(WomFloatType)
}
