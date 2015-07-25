package cromwell.binding.types

import cromwell.binding.values.{WdlFile, WdlString, WdlArray}
import spray.json.JsArray

case class WdlArrayType(memberType: WdlType) extends WdlType {
  val toWdlString: String = s"Array[${memberType.toWdlString}]"

  private def coerceIterable(s: Seq[Any]): WdlArray = {
    val coerced = s.map {memberType.coerceRawValue(_).get}
    WdlArray(WdlArrayType(coerced.head.wdlType), coerced)
  }

  override protected def coercion = {
    case s: Seq[Any] if s.nonEmpty => coerceIterable(s)
    case js: JsArray if js.elements.nonEmpty => coerceIterable(js.elements)
    case wdlArray: WdlArray => wdlArray.wdlType.memberType match {
      case WdlStringType if memberType == WdlFileType =>
        // Coerce Array[String] -> Array[File]
        WdlArray(WdlArrayType(WdlFileType), wdlArray.value.map(str => WdlFile(str.asInstanceOf[WdlString].value)).toList)
    }
  }
}
