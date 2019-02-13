package wom.types

import cats.instances.list._
import cats.syntax.traverse._
import common.util.TryUtil
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import spray.json.JsObject
import wom.values._

import scala.util.{Failure, Success, Try}

trait WomObjectTypeLike extends WomType {
  /**
    * Validate the values against the WomObjectTypeLike and return a valid map of values if possible.
    * This is an indirection from the usual coercion path but it allows WomObject to validate values against the WomObjectTypeLike at
    * instantiation time and ensure the WomObject is only built if the values match the constraints of the type.
    * See apply method in WomObject.
   */
  def validateAndCoerceValues(values: Map[String, Any]): ErrorOr[Map[String, WomValue]] = {
    values.toList.traverse({
      case (k, v) => WomAnyType.coerceRawValue(v).toErrorOr.map(k -> _)
    }).map(_.toMap)
  }
}

case object WomObjectType extends WomObjectTypeLike {
  val stableName: String = "Object"

  private def handleCoercionFailures(tries: Try[_]*) = {
    val errorMessages = tries collect {
      case Failure(f) => f.getMessage
    } mkString ","

    throw new UnsupportedOperationException(s"Coercion failed: $errorMessages")
  }

  override protected def coercion = {
    case o: WomObjectLike => WomObject(o.values)
    case m: WomMap if isMapCoercable(m) =>
      val coercedMap = m.value map {
        case (k, v) => toWomString(k) -> toWomString(v)
      } collect {
        case (Success(k), Success(v)) => k.value -> v
        case (k, v) => handleCoercionFailures(k, v)
      }

      WomObject(coercedMap)
    case js: JsObject =>

      val mapToTry = js.fields map { case (key, value) => key -> WomAnyType.coerceRawValue(value) }
      val mapOfTry = mapToTry map { kvp => kvp._2 map { kvp._1 -> _ } }
      // The TryUtil exception is ignored, we only use it to tell whether it worked or not. We use handleCoercionFailures
      // to compose the errors.
      TryUtil.sequence(mapOfTry.toList) match {
        case Success(map) => WomObject(map.toMap)
        case Failure(_) => handleCoercionFailures(mapOfTry.toSeq: _*)
      }
  }

  private def toWomString(v: WomValue) = WomStringType.coerceRawValue(v).map(_.asInstanceOf[WomString])

  override def typeSpecificIsCoerceableFrom(otherType: WomType) = otherType match {
    case _: WomObjectTypeLike => true
    case t: WomMapType if isMapTypeCoercable(t) => true
    case _ => false
  }

  def isMapTypeCoercable(t: WomMapType) = WomStringType.isCoerceableFrom(t.keyType) && WomStringType.isCoerceableFrom(t.valueType)
  def isMapCoercable(m: WomMap) = isMapTypeCoercable(m.womType)
}
