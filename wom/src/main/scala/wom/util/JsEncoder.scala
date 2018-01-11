package wom.util

import cats.data.Validated
import cats.instances.list._
import cats.syntax.apply._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import wom.types.WomNothingType
import wom.values.{WomArray, WomBoolean, WomFloat, WomInteger, WomMap, WomOptionalValue, WomString, WomValue}

import scala.collection.JavaConverters._

/**
  * Converts a WomValue into a javascript compatible value.
  */
class JsEncoder {

  /**
    * Base implementation converts any WomPrimitive (except WomFile) into a javascript compatible value.
    *
    * Inputs, and returned output must be one of:
    * - WomString
    * - WomBoolean
    * - WomFloat
    * - WomInteger
    * - WomMap
    * - WomArray
    * - A "WomNull" equal to WomOptionalValue(WomNothingType, None)
    *
    * The WomMap keys and values, and WomArray elements must be the one of the above, recursively.
    *
    * Instances of WomFile are not permitted, and must be already converted to one of the above types.
    *
    * @param value A WOM value.
    * @return The javascript equivalent.
    */
  def encode(value: WomValue): ErrorOr[AnyRef] = {
    value match {
      case WomOptionalValue(WomNothingType, None) => Validated.valid(null)
      case WomString(string) => string.valid
      case WomInteger(int) => Int.box(int).valid
      case WomFloat(double) => Double.box(double).valid
      case WomBoolean(boolean) => Boolean.box(boolean).valid
      case WomArray(_, array) => array.toList.traverse[ErrorOr, AnyRef](encode).map(JsArray)
      case WomMap(_, map) => map.traverse({
        case (mapKey, mapValue) => (encodeString(mapKey), encode(mapValue)).mapN((_, _))
      }).map(JsMap)
      case _ => s"$getClass is unable to encode value: $value".invalidNel
    }
  }

  def encodeString(value: WomValue): ErrorOr[String] = {
    encode(value) flatMap {
      case string: String => string.valid
      case other =>
        // http://2ality.com/2012/03/converting-to-string.html
        val jsString = JsUtil.evalRaw(""""" + other""", Map("other" -> other).asJava)
        jsString flatMap {
          case string: String => string.valid
          case unexpected => s"Expected to convert '$value' to a String but ended up with '$unexpected'".invalidNel
        }
    }
  }
}
