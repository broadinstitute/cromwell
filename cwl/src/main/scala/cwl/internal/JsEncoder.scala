package cwl.internal

import mouse.all._
import JsUtil.{Js, JsArray, JsObject, JsPrimitive}
import cats.data.Validated.Valid
import common.validation.ErrorOr.ErrorOr
import wom.values.{WomArray, WomBoolean, WomFloat, WomInteger, WomMap, WomObjectLike, WomOptionalValue, WomString, WomValue}

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
  def encode(value: WomValue): Js = {
    value match {
      case WomOptionalValue(_, None) => JsPrimitive(null)
      case WomOptionalValue(_, Some(innerValue)) => encode(innerValue)
      case WomString(string) => string |> JsPrimitive
      case WomInteger(int) => Int.box(int) |> JsPrimitive
      case WomFloat(double) => Double.box(double) |> JsPrimitive
      case WomBoolean(boolean) => Boolean.box(boolean) |> JsPrimitive
      case WomArray(_, array) => array.toList.map(encode).toArray |> JsArray
      case WomMap(_, map) => map.map{
        case (mapKey, mapValue) => (encodeString(mapKey), encode(mapValue))
      } |> JsObject
      case objectLike: WomObjectLike => objectLike.values.map{
        case (key, innerValue) => (key, encode(innerValue))
      } |> JsObject
      case _ => throw new RuntimeException(s"$getClass is unable to encode value: $value")
    }
  }

  def encodeString(value: WomValue): String = {
    encode(value) match {
      case JsPrimitive(string: String) => string
      case _ =>
        val jsString: ErrorOr[WomValue] = JsUtil.evalStructish(""""" + other""","other" -> value)
        jsString match {
          case Valid(WomString(string)) => string
          case unexpected => throw new RuntimeException(s"Expected to convert '$value' to a String but ended up with '$unexpected'")
        }
      /*
    case JsField(other) =>
      TODO
      // http://2ality.com/2012/03/converting-to-string.html
      */
    }
  }
}
