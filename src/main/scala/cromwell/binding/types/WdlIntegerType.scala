package cromwell.binding.types

import cromwell.binding.values.WdlInteger

case object WdlIntegerType extends WdlType {
  override def toWdlString: String = "Int"

  override protected def coercion = {
    case i: Integer => WdlInteger(i)
  }
}