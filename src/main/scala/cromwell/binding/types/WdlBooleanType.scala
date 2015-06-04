package cromwell.binding.types

import cromwell.binding.values.WdlBoolean
import spray.json.{JsBoolean, JsString}

case object WdlBooleanType extends WdlType {
  override def toWdlString: String = "Boolean"

  override protected def coercion = {
    case b: Boolean => WdlBoolean(b)
    case s: String if s.equalsIgnoreCase("true") => WdlBoolean.True
    case s: String if s.equalsIgnoreCase("false") => WdlBoolean.False
    case s: JsString if s.value.equalsIgnoreCase("true") => WdlBoolean.True
    case s: JsString if s.value.equalsIgnoreCase("false") => WdlBoolean.False
    case s: JsBoolean => WdlBoolean(s.value)
  }
}

