package cromwell.binding.types

import cromwell.binding.values.WdlInteger
import spray.json.JsNumber

case object WdlIntegerType extends WdlPrimitiveType {
  val toWdlString: String = "Int"

  override protected def coercion = {
    case i: Integer => WdlInteger(i)
    case n: JsNumber => WdlInteger(n.value.intValue())
  }

  override def fromWdlString(rawString: String) = WdlInteger(rawString.toInt)
}
