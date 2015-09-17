package cromwell.binding.types

import cromwell.binding.values.{WdlArray, WdlFile, WdlString}
import spray.json.JsArray

case class WdlArrayType(memberType: WdlType) extends WdlType {
  val toWdlString: String = s"Array[${memberType.toWdlString}]"

  private def coerceIterable(values: Seq[Any]): WdlArray = values match {
    case s:Seq[Any] if s.nonEmpty =>
      val coerced = s.map {memberType.coerceRawValue(_).get}
      WdlArray(WdlArrayType(coerced.head.wdlType), coerced)
    case _ => WdlArray(WdlArrayType(memberType), Seq())
  }

  override protected def coercion = {
    case s: Seq[Any] => coerceIterable(s)
    case js: JsArray => coerceIterable(js.elements)
    case wdlArray: WdlArray if wdlArray.wdlType.memberType == WdlStringType && memberType == WdlFileType =>
      WdlArray(WdlArrayType(WdlFileType), wdlArray.value.map(str => WdlFile(str.asInstanceOf[WdlString].value)).toList)
    case wdlArray: WdlArray if wdlArray.wdlType.memberType == memberType => wdlArray
    case wdlArray: WdlArray if wdlArray.wdlType.memberType == WdlAnyType => coerceIterable(wdlArray.value)
  }

  override def isCoerceableFrom(otherType: WdlType): Boolean = otherType match {
    case a: WdlArrayType => memberType.isCoerceableFrom(a.memberType)
    case _ => false
  }
}
