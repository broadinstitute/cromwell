package wdl4s.types
import wdl4s.values.{WdlOptionalValue, WdlValue}

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
  }

  override def isCoerceableFrom(otherType: WdlType): Boolean = otherType match {

    // Check for boxing coerceability:
    case wdlType: WdlType if memberType.isCoerceableFrom(wdlType) => true

    // Check inner value coerceability:
    case WdlOptionalType(otherMemberType) if memberType.isCoerceableFrom(otherMemberType) => true

    case _ => false
  }

  override def toWdlString: String = memberType.toWdlString + "?"
}
