package cromwell.binding.types

import cromwell.binding.values.{WdlValue, WdlArray, WdlFile, WdlString}
import spray.json.JsArray
import scala.util.Failure

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
    case wdlArray: WdlArray if wdlArray.wdlType.memberType.isInstanceOf[WdlArrayType] && memberType.isInstanceOf[WdlArrayType] =>
      val coercedElements = wdlArray.value.map(memberType.coerceRawValue)
      val coercedArray = coercedElements.collect({case x: Failure[_] => x}) match {
        case failures if failures.nonEmpty =>
          throw new UnsupportedOperationException(failures.map(f => s"Could not coerce $wdlArray: ${f.exception.getMessage}").mkString("\n"))
        case _ => coercedElements.map(_.get) // .get because they're all Success
      }
      WdlArray(WdlArrayType(memberType), coercedArray)

  }

  override def isCoerceableFrom(otherType: WdlType): Boolean = otherType match {
    case a: WdlArrayType => memberType.isCoerceableFrom(a.memberType)
    case _ => false
  }
}

object WdlArrayType {

  implicit class WdlArrayEnhanced(wdlType: WdlType) extends WdlType {

    override protected def coercion: PartialFunction[Any, WdlValue] = wdlType.coercion
    override def toWdlString: String = wdlType.toWdlString

    def isAnArrayOf(genericType: WdlType) = wdlType.isInstanceOf[WdlArrayType] && wdlType.asInstanceOf[WdlArrayType].memberType == genericType
  }
}
