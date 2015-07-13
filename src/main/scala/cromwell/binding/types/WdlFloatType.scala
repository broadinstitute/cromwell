package cromwell.binding.types

import cromwell.binding.values.WdlFloat
import spray.json.JsNumber

case object WdlFloatType extends WdlPrimitiveType {
  val toWdlString: String = "Float"

  override protected def coercion = {
    case d: Double => WdlFloat(d)
    case n: JsNumber => WdlFloat(n.value.doubleValue())
  }

  override def fromRawString(rawString: String) = WdlFloat(rawString.toFloat)
}
