package wdl4s.types

import lenthall.util.TryUtil
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlPair, WdlString, WdlValue}
import spray.json.JsArray
import wdl4s.values.WdlArray.WdlArrayLike

import scala.collection.JavaConverters._
import scala.util.{Failure, Success}

sealed trait WdlArrayType extends WdlType {

  val memberType: WdlType
  val guaranteedNonEmpty: Boolean

  private def coerceIterable(values: Seq[Any]): WdlArray = values match {
    case s:Seq[Any] if s.nonEmpty =>
      val coerced = s.map {memberType.coerceRawValue(_).get}
      WdlArray(this, coerced)
    case _ => WdlArray(this, Seq())
  }

  private val allowEmpty = !guaranteedNonEmpty
  override protected def coercion: PartialFunction[Any, WdlValue] = {
    case s: Seq[Any] if allowEmpty || s.nonEmpty => coerceIterable(s)
    case js: JsArray if allowEmpty || js.elements.nonEmpty => coerceIterable(js.elements)
    case javaList: java.util.List[_] if allowEmpty || !javaList.isEmpty => coerceIterable(javaList.asScala)
    case WdlArray(WdlMaybeEmptyArrayType.EmptyArrayType, _) => WdlArray(this, Seq.empty)
    case wdlArray: WdlArray if (allowEmpty || wdlArray.nonEmpty) && wdlArray.wdlType.memberType == WdlStringType && memberType == WdlFileType =>
      WdlArray(WdlArrayType(WdlFileType, guaranteedNonEmpty), wdlArray.value.map(str => WdlFile(str.asInstanceOf[WdlString].value)).toList)
    case wdlArray: WdlArray if (allowEmpty || wdlArray.nonEmpty) && wdlArray.wdlType.memberType == memberType => wdlArray
    case wdlArray: WdlArray if (allowEmpty || wdlArray.nonEmpty) && wdlArray.wdlType.memberType == WdlAnyType => coerceIterable(wdlArray.value)
    case wdlArray: WdlArray if (allowEmpty || wdlArray.nonEmpty) && wdlArray.wdlType.memberType.isInstanceOf[WdlArrayType] && memberType.isInstanceOf[WdlArrayType] =>
      TryUtil.sequence(wdlArray.value.map(memberType.coerceRawValue)) match {
        case Success(values) => WdlArray(this, values)
        case Failure(ex) => throw ex
      }
    case wdlArray: WdlArray if (allowEmpty || wdlArray.nonEmpty) && memberType.isCoerceableFrom(wdlArray.wdlType.memberType) =>
      wdlArray.map(v => memberType.coerceRawValue(v).get) // .get because isCoerceableFrom should make it safe
    case WdlArrayLike(wdlArray) if this.isCoerceableFrom(wdlArray.wdlType) => coercion.apply(wdlArray)
    case wdlValue: WdlValue if memberType.isCoerceableFrom(wdlValue.wdlType) =>
      memberType.coerceRawValue(wdlValue) match {
        case Success(coercedValue) => WdlArray(this, Seq(coercedValue))
        case Failure(ex) => throw ex
      }
  }

  override def isCoerceableFrom(otherType: WdlType): Boolean = otherType match {
    case WdlArrayType(otherMemberType) => memberType.isCoerceableFrom(otherMemberType) || otherMemberType == WdlNothingType
    case mapType: WdlMapType => matchesMapType(mapType)
    case _ => false
  }

  def asNonEmptyArrayType = WdlNonEmptyArrayType(memberType)

  private def matchesMapType(wdlMapType: WdlMapType) = (this, wdlMapType) match {
    case (WdlArrayType(WdlPairType(leftType, rightType)), WdlMapType(keyType, valueType)) => leftType.isCoerceableFrom(keyType) && rightType.isCoerceableFrom(valueType)
    case _ => false
  }

}

case class WdlMaybeEmptyArrayType(memberType: WdlType) extends WdlArrayType {
  override val toWdlString: String = s"Array[${memberType.toWdlString}]"
  override val guaranteedNonEmpty = false
}

object WdlMaybeEmptyArrayType {
  val EmptyArrayType = WdlMaybeEmptyArrayType(WdlNothingType)
}

case class WdlNonEmptyArrayType(memberType: WdlType) extends WdlArrayType {
  override val toWdlString: String = s"Array[${memberType.toWdlString}]+"
  override val guaranteedNonEmpty = true
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
