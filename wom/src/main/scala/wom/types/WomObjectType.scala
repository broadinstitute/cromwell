package wom.types

import spray.json.JsObject
import wom.values._

import scala.util.{Failure, Success, Try}

case object WomObjectType extends WomType {
  val toDisplayString: String = "Object"

  private def handleCoercionFailures(tries: Try[_]*) = {
    val errorMessages = tries collect {
      case Failure(f) => f.getMessage
    } mkString ","

    throw new UnsupportedOperationException(s"Coercion failed: $errorMessages")
  }

  override protected def coercion = {
    case o: WomObject => o
    case m: WomMap if isMapCoercable(m) =>
      val coercedMap = m.value map {
        case (k, v) => toWomString(k) -> toWomString(v)
      } collect {
        case (Success(k), Success(v)) => k.value -> v
        case (k, v) => handleCoercionFailures(k, v)
      }

      WomObject(coercedMap)
    case js: JsObject =>
      val coercedMap = WomMap.coerceMap(js.fields, WomMapType(WomStringType, WomAnyType)).value map {
        // get is safe because coerceMap above would have failed already if k was not coerceable to WomString
        case (k, v) => toWomString(k).get.value -> v
      }

      WomObject(coercedMap)
  }

  private def toWomString(v: WomValue) = WomStringType.coerceRawValue(v).map(_.asInstanceOf[WomString])

  override def isCoerceableFrom(otherType: WomType) = otherType match {
    case WomObjectType => true
    case t: WomMapType if isMapTypeCoercable(t) => true
    case _ => false
  }

  def isMapTypeCoercable(t: WomMapType) = WomStringType.isCoerceableFrom(t.keyType) && WomStringType.isCoerceableFrom(t.valueType)
  def isMapCoercable(m: WomMap) = isMapTypeCoercable(m.womType)
}
