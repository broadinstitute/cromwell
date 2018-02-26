package cwl.internal

import mouse.all._
import JsUtil.{ECMAScriptVariable, ESArray, ESObject, ESPrimitive}
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
  def encode(value: WomValue): ECMAScriptVariable = {
    value match {
      case WomOptionalValue(_, None) => ESPrimitive(null)
      case WomOptionalValue(_, Some(innerValue)) => encode(innerValue)
      case WomString(string) => string |> ESPrimitive
      case WomInteger(int) => Int.box(int) |> ESPrimitive
      case WomFloat(double) => Double.box(double) |> ESPrimitive
      case WomBoolean(boolean) => Boolean.box(boolean) |> ESPrimitive
      case WomArray(_, array) => array.toList.map(encode).toArray |> ESArray
      case WomMap(_, map) => map.map{
        case (mapKey, mapValue) => (encodeString(mapKey), encode(mapValue))
      } |> ESObject
      case objectLike: WomObjectLike => objectLike.values.map{
        case (key, innerValue) => (key, encode(innerValue))
      } |> ESObject
      case _ => throw new RuntimeException(s"$getClass is unable to encode value: $value")
    }
  }

  def encodeString(value: WomValue): String = {
    encode(value) match {
      case ESPrimitive(string: String) => string
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
