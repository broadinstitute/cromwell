package wom.types

import spray.json.{JsNumber, JsString}
import wom.values.{WomInteger, WomLong, WomString}

case object WomLongType extends WomPrimitiveType {
  val stableName: String = "Long"

  override protected def coercion = {
    case i: Long => WomLong(i)
    case i: Integer => WomLong(i.toLong)
    case n: JsNumber if n.value.isValidLong => WomLong(n.value.longValue())
    case WomInteger(i) => WomLong(i.toLong)
    case i: WomLong => i
    case s: WomString => WomLong(s.value.toLong)
    case s: String => WomLong(s.toLong)
    case s: JsString => WomLong(s.value.toLong)
  }
}
