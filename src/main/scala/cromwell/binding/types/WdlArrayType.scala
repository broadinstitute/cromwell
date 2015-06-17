package cromwell.binding.types

import cromwell.binding.values.{WdlInteger, WdlArray}
import spray.json.JsArray

case class WdlArrayType(memberType: WdlType) extends WdlType {
  val toWdlString: String = s"Array[${memberType.toWdlString}]"

  private def coerceIterable(s: Seq[Any]): WdlArray = {
    val coerced = s.map {e => memberType.coerceRawValue(e).get}
    WdlArray(WdlArrayType(coerced.head.wdlType), coerced)
  }

  override protected def coercion = {
    case s: Seq[Any] if s.nonEmpty => coerceIterable(s)
    case js: JsArray if js.elements.nonEmpty => coerceIterable(js.elements)
  }

  override def fromRawString(s: String) = ???
}
