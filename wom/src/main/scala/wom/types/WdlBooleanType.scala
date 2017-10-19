package wom.types

import spray.json.{JsBoolean, JsString}
import wom.values.WdlBoolean

import scala.util.{Success, Try}

case object WdlBooleanType extends WdlPrimitiveType {
  val toWdlString: String = "Boolean"

  override protected def coercion = {
    case b: Boolean => WdlBoolean(b)
    case s: String if s.equalsIgnoreCase("true") => WdlBoolean.True
    case s: String if s.equalsIgnoreCase("false") => WdlBoolean.False
    case s: JsString if s.value.equalsIgnoreCase("true") => WdlBoolean.True
    case s: JsString if s.value.equalsIgnoreCase("false") => WdlBoolean.False
    case s: JsBoolean => WdlBoolean(s.value)
    case b: WdlBoolean => b
  }

  private def comparisonOperator(rhs: WdlType, symbol: String): Try[WdlType] = rhs match {
    case WdlBooleanType => Success(WdlBooleanType)
    case WdlOptionalType(memberType) => comparisonOperator(memberType, symbol)
    case _ => invalid(s"$this $symbol $rhs")
  }

  override def equals(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, "==")
  override def lessThan(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, "<")
  override def greaterThan(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, ">")
  override def or(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, "||")
  override def and(rhs: WdlType): Try[WdlType] = comparisonOperator(rhs, "&&")
  override def not: Try[WdlType] = Success(WdlBooleanType)
}
