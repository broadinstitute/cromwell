package cromwell.binding.types

import cromwell.binding.values.WdlString

case object WdlStringType extends WdlType {
  override def toWdlString: String = "String"

  override protected def coercion = {
    case s: String => WdlString(s)
  }
}
