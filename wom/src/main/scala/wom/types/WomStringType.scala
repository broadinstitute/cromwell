package wom.types

import spray.json.JsString
import wom.WomExpressionException
import wom.values.{WomPrimitive, WomPrimitiveFile, WomString, WomValue}

import scala.util.{Failure, Success, Try}

case object WomStringType extends WomPrimitiveType {
  val toDisplayString: String = "String"

  override def coercion(): PartialFunction[Any, WomValue] = {
    case s: String => WomString(s)
    case s: JsString => WomString(s.value)
    case s: WomString => s
    case f: WomPrimitiveFile => WomString(f.value)
    case p: WomPrimitive => WomString(p.toWomString)
  }

  private def comparisonOperator(rhs: WomType, symbol: String): Try[WomType] = rhs match {
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

  override def equals(rhs: WomType): Try[WomType] =

    rhs match {
    case WomOptionalType(wct:WomCoproductType) =>
      wct.types.exists(_.equals(WomStringType)) match {
        case true => Success(WomBooleanType)
        case _ => Failure(new WomExpressionException(s"Type equality could not be asserted because $rhs was not found in the coproduct of ${wct.toDisplayString}"))
      }
    case other => comparisonOperator(other, "==")
  }
  override def lessThan(rhs: WomType): Try[WomType] = comparisonOperator(rhs, "<")
  override def greaterThan(rhs: WomType): Try[WomType] = comparisonOperator(rhs, ">")
}
