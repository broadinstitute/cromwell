package wom.types

import spray.json.{JsNumber, JsString}
import wom.types.WomIntegerLike._
import wom.values.{WomInteger, WomLong, WomString}

import scala.util.{Success, Try}

case object WomIntegerType extends WomPrimitiveType {
  val stableName: String = "Int"

  override protected def coercion = {
    case i: Integer => WomInteger(i)
    case n: JsNumber if n.value.isValidInt => WomInteger(n.value.intValue())
    case i: WomInteger => i
    case WomLong(i) if i.inIntRange => WomInteger(i.toInt)
    case WomLong(i) => throw new RuntimeException(
      s"Tried to convert a Long value $i into an Int but it was outside the bounds of acceptable Ints, namely ${Int.MinValue} <-> ${Int.MaxValue}")
    case s: WomString => {
        WomInteger(s.value.toInt)
    }
    case s: String =>
      val bigTry = Try(BigDecimal(s))
      if (bigTry.isSuccess)
        WomInteger(bigTry.get.intValue())
      else
        WomInteger(s.toInt)
    case s: JsString => WomInteger(s.value.toInt)
  }

  private def binaryOperator(rhs: WomType, symbol: String): Try[WomType] = rhs match {
    case WomIntegerType => Success(WomIntegerType)
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
  override def unaryPlus: Try[WomType] = Success(WomIntegerType)
  override def unaryMinus: Try[WomType] = Success(WomIntegerType)
}
