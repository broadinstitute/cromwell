package wdl.shared.model.expression

import spray.json.{JsArray, JsBoolean, JsNull, JsNumber, JsObject, JsString, JsValue}
import wom.TsvSerializable
import wom.values.{
  WomArray,
  WomBoolean,
  WomFile,
  WomFloat,
  WomInteger,
  WomMap,
  WomObjectLike,
  WomOptionalValue,
  WomPair,
  WomString,
  WomValue
}

import scala.util.Try

object ValueEvaluation {
  def serializeWomValue[A <: WomValue with TsvSerializable](functionName: String,
                                                            womValue: WomValue,
                                                            defaultIfOptionalEmpty: A
  ): Try[String] = {
    val wdlClass = defaultIfOptionalEmpty.getClass
    def castOrDefault(womValue: WomValue): A = womValue match {
      case WomOptionalValue(_, None) => defaultIfOptionalEmpty
      case WomOptionalValue(_, Some(v)) => wdlClass.cast(v)
      case _ => wdlClass.cast(womValue)
    }

    for {
      downcast <- Try(castOrDefault(womValue))
      tsvSerialized <- downcast.tsvSerialize
    } yield tsvSerialized
  }

  def valueToJson(womValue: WomValue): JsValue = womValue match {
    case WomInteger(i) => JsNumber(i)
    case WomFloat(f) => JsNumber(f)
    case WomString(s) => JsString(s)
    case WomBoolean(b) => JsBoolean(b)
    case f: WomFile => JsString(f.value)
    case WomPair(left, right) => JsObject(Map("left" -> valueToJson(left), "right" -> valueToJson(right)))
    case WomArray(_, values) => JsArray(values.map(valueToJson).toVector)
    case WomMap(_, value) => JsObject(value map { case (k, v) => k.valueString -> valueToJson(v) })
    case o: WomObjectLike => JsObject(o.values map { case (k, v) => k -> valueToJson(v) })
    case opt: WomOptionalValue =>
      opt.value match {
        case Some(inner) => valueToJson(inner)
        case None => JsNull
      }
  }
}
