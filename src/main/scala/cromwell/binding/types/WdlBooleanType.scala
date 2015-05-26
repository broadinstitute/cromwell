package cromwell.binding.types

import cromwell.binding.values.WdlBoolean

case object WdlBooleanType extends WdlType {
  override def toWdlString: String = "Boolean"

  override protected def coercion = {
    case b: Boolean => WdlBoolean(b)
    case s: String if s.equalsIgnoreCase("true") => WdlBoolean.True
    case s: String if s.equalsIgnoreCase("false") => WdlBoolean.False
  }
}

