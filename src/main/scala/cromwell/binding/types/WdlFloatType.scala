package cromwell.binding.types

import cromwell.binding.values.WdlFloat
import spray.json.JsNumber

case object WdlFloatType extends WdlType {
  override def toWdlString: String = "Float"

  override protected def coercion = {
    case d: Double => WdlFloat(d)
    case n: JsNumber => WdlFloat(n.value.doubleValue())
  }
}

