package cromwell.binding.types

import cromwell.binding.values.{WdlInteger, WdlString}
import spray.json.{JsNumber, JsString}

case object WdlIntegerType extends WdlPrimitiveType {
  val toWdlString: String = "Int"

  override protected def coercion = {
    case i: Integer => WdlInteger(i)
    case n: JsNumber => WdlInteger(n.value.intValue())
    case i: WdlInteger => i
    case s: WdlString => WdlInteger(s.value.toInt)
    case s: String => WdlInteger(s.toInt)
    case s: JsString => WdlInteger(s.value.toInt)
  }

  override def fromWdlString(rawString: String) = WdlInteger(rawString.toInt)
}
