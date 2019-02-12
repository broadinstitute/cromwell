package wom.types

import cats.data.NonEmptyList
import spray.json.JsString
import wom.values.{WomEnumerationValue, WomValue}


/**
  * An enumeration of possible states a value can inhabit.
  */
case class WomEnumerationType(values: NonEmptyList[String]) extends WomPrimitiveType {

  override def coercion: PartialFunction[Any, WomValue] = {
    case womValue: WomValue if (values.toList.contains(womValue.valueString)) => WomEnumerationValue(this, womValue.valueString)
    case name: String if (values.toList.contains(name)) => WomEnumerationValue(this, name)
    case JsString(name)  if (values.toList.contains(name)) => WomEnumerationValue(this, name)
  }

  override def stableName: String =
    values.toList.mkString("Enumeration[",", ", "]")
}
