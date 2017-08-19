package wdl4s.wdl.types
import spray.json.JsNull
import wdl4s.wdl.values.{WdlOptionalValue, WdlValue}

import scala.util.Try

case class WdlOptionalType(memberType: WdlType) extends WdlType {
  /**
    * Method to be overridden by implementation classes defining a partial function
    * for the conversion of raw input values to specific implementation class value types.
    * i.e.  `WdlBooleanType` should define a partial function that knows how to
    * construct `WdlBoolean`s for inputs of supported types and contents.  Values for which
    * the partial function is not defined are assumed to not be convertible to the target type.
    */
  override protected def coercion: PartialFunction[Any, WdlValue] = {

    // It's safe to box up values implicitly:
    case wdlValue: WdlValue if memberType.equals(wdlValue.wdlType) => WdlOptionalValue(wdlValue)
    case coerceable: Any if memberType.coercionDefined(coerceable) => WdlOptionalValue(memberType.coerceRawValue(coerceable).get)

    // Coercing inner values:
    case WdlOptionalValue(otherMemberType, Some(value)) if memberType.isCoerceableFrom(otherMemberType) => WdlOptionalValue(memberType, Option(memberType.coerceRawValue(value).get))
    case WdlOptionalValue(otherMemberType, None) if memberType.isCoerceableFrom(otherMemberType) => WdlOptionalValue(memberType, None)
      
    // Javascript null coerces to empty value
    case JsNull => WdlOptionalValue(memberType, None)
  }

  override def isCoerceableFrom(otherType: WdlType): Boolean = otherType match {

    // Check for boxing coerceability:
    case wdlType: WdlType if memberType.isCoerceableFrom(wdlType) => true

    // Check inner value coerceability:
    case WdlOptionalType(otherMemberType) if memberType.isCoerceableFrom(otherMemberType) => true

    case _ => false
  }

  override def add(rhs: WdlType): Try[WdlType] = memberType.add(rhs)
  override def subtract(rhs: WdlType): Try[WdlType] = memberType.subtract(rhs)
  override def multiply(rhs: WdlType): Try[WdlType] = memberType.multiply(rhs)
  override def divide(rhs: WdlType): Try[WdlType] = memberType.divide(rhs)
  override def mod(rhs: WdlType): Try[WdlType] = memberType.mod(rhs)
  override def equals(rhs: WdlType): Try[WdlType] = memberType.equals(rhs)
  override def lessThan(rhs: WdlType): Try[WdlType] = memberType.lessThan(rhs)
  override def greaterThan(rhs: WdlType): Try[WdlType] = memberType.greaterThan(rhs)
  override def or(rhs: WdlType): Try[WdlType] = memberType.or(rhs)
  override def and(rhs: WdlType): Try[WdlType] = memberType.and(rhs)
  override def not: Try[WdlType] = memberType.not
  override def unaryPlus: Try[WdlType] = memberType.unaryPlus
  override def unaryMinus: Try[WdlType] = memberType.unaryMinus

  override def toWdlString: String = memberType.toWdlString + "?"

  def none = WdlOptionalValue.none(memberType)
}
