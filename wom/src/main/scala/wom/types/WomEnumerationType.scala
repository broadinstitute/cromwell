package wom.types

import cats.data.NonEmptyList
import wom.values.{WomEnumerationValue, WomValue}


/**
  * An enumeration of possible states a value can inhabit.
  */
case class WomEnumerationType(values: NonEmptyList[String]) extends WomPrimitiveType {

  override def coercion: PartialFunction[Any, WomValue] = {
    case womValue: WomValue if (values.toList.contains(womValue.valueString)) => WomEnumerationValue(this, womValue.valueString)
  }

  override def toDisplayString: String =
    values.toList.mkString("Enumeration[",", ", "]")
}
