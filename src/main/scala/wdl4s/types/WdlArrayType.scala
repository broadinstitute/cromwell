package wdl4s.types

import lenthall.util.TryUtil
import wdl4s.values.{WdlArray, WdlFile, WdlString, WdlValue}
import spray.json.JsArray

import scala.collection.JavaConverters._
import scala.util.{Failure, Success}

sealed trait WdlArrayType extends WdlType {

  val memberType: WdlType
  val guaranteedNonEmpty: Boolean

  private def coerceIterable(values: Seq[Any]): WdlArray = values match {
    case s:Seq[Any] if s.nonEmpty =>
      val coerced = s.map {memberType.coerceRawValue(_).get}
      WdlArray(WdlArrayType(coerced.head.wdlType), coerced)
    case _ => WdlArray(WdlArrayType(memberType), Seq())
  }

  protected def coercionMaker(allowEmpty: Boolean): PartialFunction[Any, WdlValue] = {
    case s: Seq[Any] if allowEmpty || s.nonEmpty => coerceIterable(s)
    case js: JsArray if allowEmpty || js.elements.nonEmpty => coerceIterable(js.elements)
    case javaList: java.util.List[_] if allowEmpty || !javaList.isEmpty => coerceIterable(javaList.asScala)
    case WdlArray(WdlMaybeEmptyArrayType.EmptyArrayType, _) => WdlArray(WdlArrayType(memberType), Seq.empty)
    case wdlArray: WdlArray if (allowEmpty || wdlArray.nonEmpty) && wdlArray.wdlType.memberType == WdlStringType && memberType == WdlFileType =>
      WdlArray(WdlArrayType(WdlFileType), wdlArray.value.map(str => WdlFile(str.asInstanceOf[WdlString].value)).toList)
    case wdlArray: WdlArray if (allowEmpty || wdlArray.nonEmpty) && wdlArray.wdlType.memberType == memberType => wdlArray
    case wdlArray: WdlArray if (allowEmpty || wdlArray.nonEmpty) && wdlArray.wdlType.memberType == WdlAnyType => coerceIterable(wdlArray.value)
    case wdlArray: WdlArray if (allowEmpty || wdlArray.nonEmpty) && wdlArray.wdlType.memberType.isInstanceOf[WdlArrayType] && memberType.isInstanceOf[WdlArrayType] =>
      TryUtil.sequence(wdlArray.value.map(memberType.coerceRawValue)) match {
        case Success(values) => WdlArray(WdlArrayType(memberType), values)
        case Failure(ex) => throw ex
      }
    case wdlArray: WdlArray if (allowEmpty || wdlArray.nonEmpty) && memberType.isCoerceableFrom(wdlArray.wdlType.memberType) =>
      wdlArray.map(v => memberType.coerceRawValue(v).get) // .get because isCoerceableFrom should make it safe
    case wdlValue: WdlValue if memberType.isCoerceableFrom(wdlValue.wdlType) =>
      memberType.coerceRawValue(wdlValue) match {
        case Success(coercedValue) => WdlArray(this, Seq(coercedValue))
        case Failure(ex) => throw ex
      }
  }

  override def isCoerceableFrom(otherType: WdlType): Boolean = otherType match {
    case WdlArrayType(otherMemberType) => memberType.isCoerceableFrom(otherMemberType) || otherMemberType == WdlNothingType
    case _ => false
  }

  def asNonEmptyArrayType = WdlNonEmptyArrayType(memberType)
}

case class WdlMaybeEmptyArrayType(memberType: WdlType) extends WdlArrayType {
  override val toWdlString: String = s"Array[${memberType.toWdlString}]"
  override val guaranteedNonEmpty = false
  override protected def coercion = coercionMaker(allowEmpty = true)
}

object WdlMaybeEmptyArrayType {
  val EmptyArrayType = WdlMaybeEmptyArrayType(WdlNothingType)
}

case class WdlNonEmptyArrayType(memberType: WdlType) extends WdlArrayType {
  override val toWdlString: String = s"Array[${memberType.toWdlString}]+"
  override val guaranteedNonEmpty = true
  override protected def coercion = coercionMaker(allowEmpty = false)
}

object WdlArrayType {

  def apply(memberType: WdlType, guaranteedNonEmpty: Boolean = false): WdlArrayType =
    if (guaranteedNonEmpty) WdlNonEmptyArrayType(memberType)
    else WdlMaybeEmptyArrayType(memberType)
  def unapply(wdlType: WdlType): Option[WdlType] = {
    wdlType match {
      case arrayType: WdlArrayType => Option(arrayType.memberType)
      case _ => None
    }
  }
}
