package cromwell.binding.types

import cromwell.binding.values.WdlUri
import spray.json.JsString

// FIXME: This is a placeholder until I get Chris' code
case object WdlUriType extends WdlType {
  val toWdlString: String = "URI"

  override def fromRawString(rawString: String) = WdlUri(rawString)

  override protected def coercion = {
    case s: String => WdlUri(s)
    case s: JsString => WdlUri(s.value)
  }
}
