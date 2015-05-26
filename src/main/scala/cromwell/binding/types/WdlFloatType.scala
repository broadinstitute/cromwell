package cromwell.binding.types

import cromwell.binding.values.WdlFloat

case object WdlFloatType extends WdlType {
  override def toWdlString: String = "Float"

  override protected def coercion = {
    case d: Double => WdlFloat(d)
  }
}

