package wom.types

import spray.json.JsString
import wom.values.{WomPrimitive, WomPrimitiveFile, WomString}

import scala.util.{Success, Try}

case object WomStringType extends WomPrimitiveType {
  val stableName: String = "String"

  override protected def coercion = {
    case s: String => WomString(s)
    case s: JsString => WomString(s.value)
    case s: WomString => s
    case f: WomPrimitiveFile => WomString(f.value)
    case p: WomPrimitive => WomString(p.toWomString)
  }

  private def comparisonOperator(rhs: WomType, symbol: String): Try[WomType] = rhs match {
    case wct:WomCoproductType => wct.typeExists(WomStringType)
    case WomStringType => Success(WomBooleanType)
    case WomOptionalType(memberType) => comparisonOperator(memberType, symbol)
    case _ => invalid(s"$this $symbol $rhs")
  }

  override def add(rhs: WomType): Try[WomType] = rhs match {
    case WomStringType | WomIntegerType | WomFloatType => Success(WomStringType)
    case _: WomPrimitiveFileType => Success(WomStringType)
    case WomOptionalType(memberType) => add(memberType)
    case _ => invalid(s"$this + $rhs")
  }

  override def equalsType(rhs: WomType): Try[WomType] = comparisonOperator(rhs, "==")
  override def lessThan(rhs: WomType): Try[WomType] = comparisonOperator(rhs, "<")
  override def greaterThan(rhs: WomType): Try[WomType] = comparisonOperator(rhs, ">")
}
