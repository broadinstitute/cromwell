package wdl4s.wdl.types

import wdl4s.wdl.values.{WdlFloat, WdlString}
import spray.json.{JsNumber, JsString}

import scala.util.{Try, Success}

case object WdlFloatType extends WdlPrimitiveType {
  val toWdlString: String = "Float"

  override protected def coercion = {
    case d: Double => WdlFloat(d)
    case n: JsNumber => WdlFloat(n.value.doubleValue())
    case f: WdlFloat => f
    case s: String => WdlFloat(s.toDouble)
    case s: JsString => WdlFloat(s.value.toDouble)
    case s: WdlString => WdlFloat(s.value.toDouble)
  }

  override def fromWdlString(rawString: String) = WdlFloat(rawString.toDouble)

  private def binaryOperator(rhs: WdlType, symbol: String): Try[WdlType] = rhs match {
    case WdlIntegerType => Success(WdlFloatType)
    case WdlFloatType => Success(WdlFloatType)
    case WdlOptionalType(memberType) => binaryOperator(memberType, symbol)
    case _ => invalid(s"$this $symbol $rhs")
  }

  private def comparisonOperator(rhs: WdlType, symbol: String): Try[WdlType] = rhs match {
    case WdlIntegerType => Success(WdlBooleanType)
    case WdlFloatType => Success(WdlBooleanType)
    case WdlOptionalType(memberType) => comparisonOperator(memberType, symbol)
    case _ => invalid(s"$this $symbol $rhs")
  }

  override def add(rhs: WdlType): Try[WdlType] = rhs match {
    case WdlStringType => Success(WdlStringType)
    case WdlOptionalType(memberType) => add(memberType)
    case t => binaryOperator(t, "+")
  }

  override def subtract(rhs: WdlType): Try[WdlType] = binaryOperator(rhs, "-")
  override def multiply(rhs: WdlType): Try[WdlType] = binaryOperator(rhs, "*")
  override def divide(rhs: WdlType): Try[WdlType] = binaryOperator(rhs, "/")
  override def mod(rhs: WdlType): Try[WdlType] = binaryOperator(rhs, "%")
  override def equals(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, "==")
  override def lessThan(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, "<")
  override def greaterThan(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, ">")
  override def unaryPlus: Try[WdlType] = Success(WdlFloatType)
  override def unaryMinus: Try[WdlType] = Success(WdlFloatType)
}
