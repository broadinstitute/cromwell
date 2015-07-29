package cromwell.binding.types

import cromwell.binding.values.WdlString
import spray.json.JsString

case object WdlStringType extends WdlPrimitiveType {
  val toWdlString: String = "String"

  override protected def coercion = {
    case s: String => WdlString(s)
    case s: JsString => WdlString(s.value)
    case s: WdlString => s
  }
}
