package wom.types

import common.util.TryUtil
import spray.json.JsArray
import wom.values.WomArray.WomArrayLike
import wom.values.{WomArray, WomSingleFile, WomString, WomValue}

import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}

sealed trait WomArrayType extends WomType {

  val memberType: WomType
  val guaranteedNonEmpty: Boolean

  private def coerceIterable(values: Seq[Any]): WomArray = values match {
    case s:Seq[Any] if s.nonEmpty =>
      val coerced = s.map {memberType.coerceRawValue(_).get}
      WomArray(this, coerced)
    case _ => WomArray(this, Seq())
  }

  private val allowEmpty = !guaranteedNonEmpty
  override protected def coercion: PartialFunction[Any, WomValue] = {
    case s: Seq[Any] if allowEmpty || s.nonEmpty => coerceIterable(s)
    case js: JsArray if allowEmpty || js.elements.nonEmpty => coerceIterable(js.elements)
    case javaList: java.util.List[_] if allowEmpty || !javaList.isEmpty => coerceIterable(javaList.asScala)
    case WomArray(WomMaybeEmptyArrayType.EmptyArrayType, _) => WomArray(this, Seq.empty)
    case womArray: WomArray
      if (allowEmpty || womArray.nonEmpty)
        && womArray.womType.memberType == WomStringType
        && memberType == WomSingleFileType =>
      WomArray(this, womArray.value.map(str => WomSingleFile(str.asInstanceOf[WomString].value)).toList)
    case womArray: WomArray if (allowEmpty || womArray.nonEmpty) && womArray.womType.memberType == memberType => WomArray(this, womArray.value)
    case womArray: WomArray if (allowEmpty || womArray.nonEmpty) && womArray.womType.memberType == WomAnyType => coerceIterable(womArray.value)
    case womArray: WomArray if (allowEmpty || womArray.nonEmpty) && womArray.womType.memberType.isInstanceOf[WomArrayType] && memberType.isInstanceOf[WomArrayType] =>
      TryUtil.sequence(womArray.value.map(memberType.coerceRawValue)) match {
        case Success(values) => WomArray(this, values)
        case Failure(ex) => throw ex
      }
    case womArray: WomArray if (allowEmpty || womArray.nonEmpty) && memberType.isCoerceableFrom(womArray.womType.memberType) =>
      womArray.map(v => memberType.coerceRawValue(v).get) // .get because isCoerceableFrom should make it safe
    case WomArrayLike(womArray) if this.isCoerceableFrom(womArray.womType) => coercion.apply(womArray)
  }

  override def typeSpecificIsCoerceableFrom(otherType: WomType): Boolean = otherType match {
    case WomArrayType(otherMemberType) => memberType.isCoerceableFrom(otherMemberType)
    case mapType: WomMapType => isCoerceableFrom(mapType.equivalentArrayType)
    case _ => false
  }

  def asNonEmptyArrayType = WomNonEmptyArrayType(memberType)
}

case class WomMaybeEmptyArrayType(memberType: WomType) extends WomArrayType {
  override val stableName: String = s"Array[${memberType.stableName}]"
  override val guaranteedNonEmpty = false

  override def equalsType(rhs: WomType): Try[WomType] = rhs match {
    case WomMaybeEmptyArrayType(`memberType`) => Success(WomBooleanType)
    case _ => Failure(new RuntimeException(s"type $rhs was incompatible with $this"))
  }
}

object WomMaybeEmptyArrayType {
  val EmptyArrayType = WomMaybeEmptyArrayType(WomNothingType)
}

case class WomNonEmptyArrayType(memberType: WomType) extends WomArrayType {
  override val stableName: String = s"Array[${memberType.stableName}]+"
  override val guaranteedNonEmpty = true

  override def equalsType(rhs: WomType): Try[WomType] = rhs match {
    case WomNonEmptyArrayType(`memberType`) => Success(WomBooleanType)
    case _ => Failure(new RuntimeException(s"type $rhs was incompatible with $this"))
  }
}

object WomArrayType {
  def apply(memberType: WomType, guaranteedNonEmpty: Boolean = false): WomArrayType =
    if (guaranteedNonEmpty) WomNonEmptyArrayType(memberType)
    else WomMaybeEmptyArrayType(memberType)

  def unapply(womType: WomType): Option[WomType] = {
    womType match {
      case arrayType: WomArrayType => Option(arrayType.memberType)
      case mapType: WomMapType => Option(mapType.equivalentArrayType.memberType)
      case _ => None
    }
  }
}
