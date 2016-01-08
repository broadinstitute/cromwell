package wdl4s.types

import wdl4s.values.{WdlInteger, WdlString}
import spray.json.{JsNumber, JsString}

import scala.util.{Try, Success}

case object WdlIntegerType extends WdlPrimitiveType {
  val toWdlString: String = "Int"

  override protected def coercion = {
    case i: Integer => WdlInteger(i)
    case n: JsNumber => WdlInteger(n.value.intValue())
    case i: WdlInteger => i
    case s: WdlString => WdlInteger(s.value.toInt)
    case s: String => WdlInteger(s.toInt)
    case s: JsString => WdlInteger(s.value.toInt)
  }

  override def fromWdlString(rawString: String) = WdlInteger(rawString.toInt)

  private def binaryOperator(rhs: WdlType, symbol: String): Try[WdlType] = rhs match {
    case WdlIntegerType => Success(WdlIntegerType)
    case WdlFloatType => Success(WdlFloatType)
    case _ => invalid(s"$this $symbol $rhs")
  }

  private def comparisonOperator(rhs: WdlType, symbol: String): Try[WdlType] = rhs match {
    case WdlIntegerType => Success(WdlBooleanType)
    case WdlFloatType => Success(WdlBooleanType)
    case _ => invalid(s"$this $symbol $rhs")
  }


  override def add(rhs: WdlType): Try[WdlType] = rhs match {
    case WdlStringType => Success(WdlStringType)
    case t => binaryOperator(t, "+")
  }

  override def subtract(rhs: WdlType): Try[WdlType] = binaryOperator(rhs, "-")
  override def multiply(rhs: WdlType): Try[WdlType] = binaryOperator(rhs, "*")
  override def divide(rhs: WdlType): Try[WdlType] = binaryOperator(rhs, "/")
  override def mod(rhs: WdlType): Try[WdlType] = binaryOperator(rhs, "%")
  override def equals(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, "==")
  override def lessThan(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, "<")
  override def greaterThan(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, ">")
  override def unaryPlus: Try[WdlType] = Success(WdlIntegerType)
  override def unaryMinus: Try[WdlType] = Success(WdlIntegerType)
}
