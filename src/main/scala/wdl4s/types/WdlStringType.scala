package wdl4s.types

import wdl4s.values.{WdlFile, WdlFloat, WdlInteger, WdlPrimitive, WdlString}
import spray.json.JsString

import scala.util.{Success, Try}

case object WdlStringType extends WdlPrimitiveType {
  val toWdlString: String = "String"

  override protected def coercion = {
    case s: String => WdlString(s)
    case s: JsString => WdlString(s.value)
    case s: WdlString => s
    case f: WdlFile => WdlString(f.value)
    case p: WdlPrimitive => WdlString(p.toWdlString)
  }

  private def comparisonOperator(rhs: WdlType, symbol: String): Try[WdlType] = rhs match {
    case WdlStringType => Success(WdlBooleanType)
    case _ => invalid(s"$this $symbol $rhs")
  }

  override def add(rhs: WdlType): Try[WdlType] = rhs match {
    case WdlStringType | WdlIntegerType | WdlFloatType | WdlFileType => Success(WdlStringType)
    case _ => invalid(s"$this + $rhs")
  }

  override def equals(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, "==")
  override def lessThan(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, "<")
  override def greaterThan(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, ">")
}
