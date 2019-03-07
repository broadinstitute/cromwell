package wom.types

import spray.json.JsNull
import wom.values.{WomOptionalValue, WomValue}

import scala.util.Try

case class WomOptionalType(memberType: WomType) extends WomType {

  def depth: Int = memberType match {
    case recursive: WomOptionalType => 1 + recursive.depth
    case _ => 1
  }
  /**
    * Method to be overridden by implementation classes defining a partial function
    * for the conversion of raw input values to specific implementation class value types.
    * i.e.  `WomBooleanType` should define a partial function that knows how to
    * construct `WomBoolean`s for inputs of supported types and contents.  Values for which
    * the partial function is not defined are assumed to not be convertible to the target type.
    */
  override protected def coercion: PartialFunction[Any, WomValue] = {

    // Javascript null coerces to empty value
    case JsNull => WomOptionalValue(memberType, None)
    case None => WomOptionalValue(memberType, None)

    // Coerce and adjust nesting level of equivalent nested conditionals:
    case womOptional: WomOptionalValue if baseMemberType.isCoerceableFrom(womOptional.womType.baseMemberType) => womOptional.coerceAndSetNestingLevel(this).get

    // It's safe to box up values implicitly:
    case womValue: WomValue if baseMemberType.isCoerceableFrom(womValue.womType) => WomOptionalValue(womValue).coerceAndSetNestingLevel(this).get

    case WomOptionalValue(WomNothingType, None) => WomOptionalValue(memberType, None)

    case null => WomOptionalValue(memberType, None)

    case coerceable: Any if baseMemberType.coercionDefined(coerceable) => WomOptionalValue(baseMemberType.coerceRawValue(coerceable).get).coerceAndSetNestingLevel(this).get
  }

  override def typeSpecificIsCoerceableFrom(otherType: WomType): Boolean = otherType match {

    // Check for boxing coerceability:
    case womType: WomType if memberType.isCoerceableFrom(womType) => true

    // Check inner value coerceability:
    case WomOptionalType(otherMemberType) if memberType.isCoerceableFrom(otherMemberType) => true

    // Check flattening:
    case WomOptionalType(otherMemberType: WomOptionalType) => baseMemberType.isCoerceableFrom(otherMemberType.baseMemberType)

    case _ => false
  }

  /**
    * Unpack any number of layers of optional-ness to get to the base member type:
    */
  def baseMemberType: WomType = flatOptionalType.memberType

  /**
    * Flatten this optional type (eg Int?[...]? to Int?)
    */
  def flatOptionalType: WomOptionalType = memberType match {
    case innerOptionalType: WomOptionalType => innerOptionalType.flatOptionalType
    case _ => this
  }

  def baseMemberTypeIsCompatibleWith(otherType: WomType): Boolean =
    baseMemberType.equalsType(otherType).isSuccess

  override def add(rhs: WomType): Try[WomType] = memberType.add(rhs)
  override def subtract(rhs: WomType): Try[WomType] = memberType.subtract(rhs)
  override def multiply(rhs: WomType): Try[WomType] = memberType.multiply(rhs)
  override def divide(rhs: WomType): Try[WomType] = memberType.divide(rhs)
  override def mod(rhs: WomType): Try[WomType] = memberType.mod(rhs)
  override def equalsType(rhs: WomType): Try[WomType] = memberType.equalsType(rhs)
  override def lessThan(rhs: WomType): Try[WomType] = memberType.lessThan(rhs)
  override def greaterThan(rhs: WomType): Try[WomType] = memberType.greaterThan(rhs)
  override def or(rhs: WomType): Try[WomType] = memberType.or(rhs)
  override def and(rhs: WomType): Try[WomType] = memberType.and(rhs)
  override def not: Try[WomType] = memberType.not
  override def unaryPlus: Try[WomType] = memberType.unaryPlus
  override def unaryMinus: Try[WomType] = memberType.unaryMinus

  override def stableName: String = memberType.stableName + "?"

  def none = WomOptionalValue.none(memberType)
}
