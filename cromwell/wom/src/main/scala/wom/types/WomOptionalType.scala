package wom.types

import spray.json.JsNull
import wom.values.{WomOptionalValue, WomValue}

import scala.util.Try

case class WomOptionalType(memberType: WomType) extends WomType {
  /**
    * Method to be overridden by implementation classes defining a partial function
    * for the conversion of raw input values to specific implementation class value types.
    * i.e.  `WomBooleanType` should define a partial function that knows how to
    * construct `WomBoolean`s for inputs of supported types and contents.  Values for which
    * the partial function is not defined are assumed to not be convertible to the target type.
    */
  override protected def coercion: PartialFunction[Any, WomValue] = {

    // It's safe to box up values implicitly:
    case womValue: WomValue if memberType.equals(womValue.womType) => WomOptionalValue(womValue)
    case coerceable: Any if memberType.coercionDefined(coerceable) => WomOptionalValue(memberType.coerceRawValue(coerceable).get)

    // Coercing inner values:
    case WomOptionalValue(otherMemberType, Some(value)) if memberType.isCoerceableFrom(otherMemberType) => WomOptionalValue(memberType, Option(memberType.coerceRawValue(value).get))
    case WomOptionalValue(otherMemberType, None) if memberType.isCoerceableFrom(otherMemberType) => WomOptionalValue(memberType, None)
      
    // Javascript null coerces to empty value
    case JsNull => WomOptionalValue(memberType, None)
  }

  override def isCoerceableFrom(otherType: WomType): Boolean = otherType match {

    // Check for boxing coerceability:
    case womType: WomType if memberType.isCoerceableFrom(womType) => true

    // Check inner value coerceability:
    case WomOptionalType(otherMemberType) if memberType.isCoerceableFrom(otherMemberType) => true

    case _ => false
  }

  override def add(rhs: WomType): Try[WomType] = memberType.add(rhs)
  override def subtract(rhs: WomType): Try[WomType] = memberType.subtract(rhs)
  override def multiply(rhs: WomType): Try[WomType] = memberType.multiply(rhs)
  override def divide(rhs: WomType): Try[WomType] = memberType.divide(rhs)
  override def mod(rhs: WomType): Try[WomType] = memberType.mod(rhs)
  override def equals(rhs: WomType): Try[WomType] = memberType.equals(rhs)
  override def lessThan(rhs: WomType): Try[WomType] = memberType.lessThan(rhs)
  override def greaterThan(rhs: WomType): Try[WomType] = memberType.greaterThan(rhs)
  override def or(rhs: WomType): Try[WomType] = memberType.or(rhs)
  override def and(rhs: WomType): Try[WomType] = memberType.and(rhs)
  override def not: Try[WomType] = memberType.not
  override def unaryPlus: Try[WomType] = memberType.unaryPlus
  override def unaryMinus: Try[WomType] = memberType.unaryMinus

  override def toDisplayString: String = memberType.toDisplayString + "?"

  def none = WomOptionalValue.none(memberType)
}
