package wom.types

import spray.json.{JsNumber, JsString}
import wom.values.{WomBigDecimal, WomInteger, WomLong, WomString}

import scala.util.{Success, Try}

case object WomBigDecimalType extends WomPrimitiveType {
  val stableName: String = "BigDecimal"

  override protected def coercion = {
    //Somewhy Integer is not implicitly conversed to BigDecimal
    case i: Integer => WomBigDecimal(BigDecimal(i))
    case bd: BigDecimal => WomBigDecimal(bd)
    case s: String => WomBigDecimal(BigDecimal(s))
    case JsNumber(n) => WomBigDecimal(n)
    case WomInteger(i) => WomBigDecimal(i)
    case WomBigDecimal(bd) => WomBigDecimal(bd)
    case WomLong(i) => WomBigDecimal(i)
    case WomString(s) => WomBigDecimal(BigDecimal(s))
    case JsString(s) => WomBigDecimal(BigDecimal(s))
  }

  private def binaryOperator(rhs: WomType, symbol: String): Try[WomType] = rhs match {
    case WomIntegerType => Success(WomBigDecimalType)
    case WomBigDecimalType => Success(WomBigDecimalType)
    case WomFloatType => Success(WomBigDecimalType)
    case WomOptionalType(_) => binaryOperator(WomBigDecimalType, symbol)
    case _ => invalid(s"$this $symbol $rhs")
  }

  private def comparisonOperator(rhs: WomType, symbol: String): Try[WomType] = rhs match {
    case wct: WomCoproductType => wct.typeExists(WomStringType)
    case WomIntegerType => Success(WomBooleanType)
    case WomBigDecimalType => Success(WomBooleanType)
    case WomFloatType => Success(WomBooleanType)
    case WomOptionalType(memberType) => comparisonOperator(memberType, symbol)
    case _ => invalid(s"$this $symbol $rhs")
  }


  override def add(rhs: WomType): Try[WomType] = rhs match {
    case WomStringType => Success(WomStringType)
    case WomOptionalType(_) => add(WomBigDecimalType)
    case t => binaryOperator(t, "+")
  }

  override def subtract(rhs: WomType): Try[WomType] = binaryOperator(rhs, "-")

  override def multiply(rhs: WomType): Try[WomType] = binaryOperator(rhs, "*")

  override def divide(rhs: WomType): Try[WomType] = binaryOperator(rhs, "/")

  override def mod(rhs: WomType): Try[WomType] = binaryOperator(rhs, "%")

  override def equalsType(rhs: WomType): Try[WomType] = comparisonOperator(rhs, "==")

  override def lessThan(rhs: WomType): Try[WomType] = comparisonOperator(rhs, "<")

  override def greaterThan(rhs: WomType): Try[WomType] = comparisonOperator(rhs, ">")

  override def unaryPlus: Try[WomType] = Success(WomBigDecimalType)

  override def unaryMinus: Try[WomType] = Success(WomBigDecimalType)
}
