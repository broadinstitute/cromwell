package wdl4s.wdl.types

import spray.json.JsString
import wdl4s.wdl.values.{WdlFile, WdlPrimitive, WdlString}

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
    case WdlOptionalType(memberType) => comparisonOperator(memberType, symbol)
    case _ => invalid(s"$this $symbol $rhs")
  }

  override def add(rhs: WdlType): Try[WdlType] = rhs match {
    case WdlStringType | WdlIntegerType | WdlFloatType | WdlFileType => Success(WdlStringType)
    case WdlOptionalType(memberType) => add(memberType)
    case _ => invalid(s"$this + $rhs")
  }

  override def equals(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, "==")
  override def lessThan(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, "<")
  override def greaterThan(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, ">")
}
