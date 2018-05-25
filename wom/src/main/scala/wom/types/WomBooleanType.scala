package wom.types

import spray.json.{JsBoolean, JsString}
import wom.values.WomBoolean

import scala.util.{Success, Try}

case object WomBooleanType extends WomPrimitiveType {
  val toDisplayString: String = "Boolean"

  override protected def coercion = {
    case b: Boolean => WomBoolean(b)
    case s: String if s.equalsIgnoreCase("true") => WomBoolean.True
    case s: String if s.equalsIgnoreCase("false") => WomBoolean.False
    case s: JsString if s.value.equalsIgnoreCase("true") => WomBoolean.True
    case s: JsString if s.value.equalsIgnoreCase("false") => WomBoolean.False
    case s: JsBoolean => WomBoolean(s.value)
    case b: WomBoolean => b
  }

  private def comparisonOperator(rhs: WomType, symbol: String): Try[WomType] = rhs match {
    case wct:WomCoproductType => wct.typeExists(WomStringType)
    case WomBooleanType => Success(WomBooleanType)
    case WomOptionalType(memberType) => comparisonOperator(memberType, symbol)
    case _ => invalid(s"$this $symbol $rhs")
  }

  override def equalsType(rhs: WomType): Try[WomType] = comparisonOperator(rhs, "==")
  override def lessThan(rhs: WomType): Try[WomType] = comparisonOperator(rhs, "<")
  override def greaterThan(rhs: WomType): Try[WomType] = comparisonOperator(rhs, ">")
  override def or(rhs: WomType): Try[WomType] = comparisonOperator(rhs, "||")
  override def and(rhs: WomType): Try[WomType] = comparisonOperator(rhs, "&&")
  override def not: Try[WomType] = Success(WomBooleanType)
}
