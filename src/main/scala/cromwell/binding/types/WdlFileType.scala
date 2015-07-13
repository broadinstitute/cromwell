package cromwell.binding.types

import cromwell.binding.values.WdlFile
import spray.json.JsString

case object WdlFileType extends WdlType {
  val toWdlString: String = "File"

  override protected def coercion = {
    case s: String => WdlFile(s)
    case s: JsString => WdlFile(s.value)
  }

  override def fromRawString(rawString: String) = WdlFile(rawString)
}
