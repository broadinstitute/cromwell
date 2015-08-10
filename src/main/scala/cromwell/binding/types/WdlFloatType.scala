package cromwell.binding.types

import cromwell.binding.values.{WdlFloat, WdlString}
import spray.json.{JsNumber, JsString}

case object WdlFloatType extends WdlPrimitiveType {
  val toWdlString: String = "Float"

  override protected def coercion = {
    case d: Double => WdlFloat(d)
    case n: JsNumber => WdlFloat(n.value.doubleValue())
    case f: WdlFloat => f
    case s: String => WdlFloat(s.toDouble)
    case s: JsString => WdlFloat(s.value.toDouble)
    case s: WdlString => WdlFloat(s.value.toDouble)
  }

  override def fromWdlString(rawString: String) = WdlFloat(rawString.toFloat)
}
