package cromwell.binding.types

import cromwell.binding.values.WdlFile

case object WdlFileType extends WdlType {
  override def toWdlString: String = "File"

  override protected def coercion = {
    case s: String => WdlFile(s)
  }
}
