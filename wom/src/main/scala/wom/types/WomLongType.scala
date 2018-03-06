package wom.types

import spray.json.{JsNumber, JsString}
import wom.values.{WomInteger, WomLong, WomString}

case object WomLongType extends WomPrimitiveType {
  val toDisplayString: String = "Long"

  override protected def coercion = {
    case i: Long => WomLong(i)
    case i: Integer => WomLong(i.toLong)
    case n: JsNumber if n.value.isValidLong => WomInteger(n.value.intValue())
    case i: WomInteger => i
    case s: WomString => WomInteger(s.value.toInt)
    case s: String => WomInteger(s.toInt)
    case s: JsString => WomInteger(s.value.toInt)
  }
}
