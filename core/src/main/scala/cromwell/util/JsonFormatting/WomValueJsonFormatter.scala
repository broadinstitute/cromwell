package cromwell.util.JsonFormatting

import spray.json._
import wom.types._
import wom.values._

object WomValueJsonFormatter extends DefaultJsonProtocol {
  implicit object WomValueJsonFormat extends RootJsonFormat[WomValue] {
    def write(value: WomValue): JsValue = value match {
      case s: WomString => JsString(s.value)
      case i: WomInteger => JsNumber(i.value)
      case f: WomFloat => JsNumber(f.value)
      case b: WomBoolean => JsBoolean(b.value)
      case f: WomSingleFile => JsString(f.value)
      case o: WomObjectLike => new JsObject(o.values map {case(k, v) => k -> write(v)})
      case a: WomArray => new JsArray(a.value.map(write).toVector)
      case m: WomMap => new JsObject(m.value map {case(k,v) => k.valueString -> write(v)})
      case q: WomPair => new JsObject(Map("left" -> write(q.left), "right" -> write(q.right)))
      case WomOptionalValue(_, Some(innerValue)) => write(innerValue)
      case WomOptionalValue(_, None) => JsNull
      case WomCoproductValue(_, innerValue) => write(innerValue)
        // handles WdlExpression
      case v: WomValue => JsString(v.toWomString)

    }

    // NOTE: This assumes a map's keys are strings. Since we're coming from JSON this is fine.
    // This won't support a map with complex keys (e.g. WomMapType(WomMapType(WomIntegerType, WomIntegerType), WomStringType)
    // That would require a more inventive JSON which splits out the key and value as full fields in their own right...
    // In addition, we make a lot of assumptions about what type of WomValue to create. Oh well... it should all fall out in the coercion (fingercrossed)!
    def read(value: JsValue): WomValue = value match {
      case JsObject(fields) =>
        val wdlFields: Map[WomValue, WomValue] = fields map {case (k, v) => WomString(k) -> read(v)}
        if (fields.isEmpty) WomMap(WomMapType(WomStringType, WomStringType), Map.empty[WomValue, WomValue])
        else WomMap(WomMapType(wdlFields.head._1.womType, wdlFields.head._2.womType), wdlFields)
      case JsArray(vector) if vector.nonEmpty => WomArray(WomArrayType(read(vector.head).womType), vector map read)
      case JsString(str) => WomString(str)
      case JsBoolean(bool) => WomBoolean(bool)
      case JsNumber(decimal) if decimal.isValidInt => WomInteger(decimal.toIntExact)
      case JsNumber(decimal) => WomFloat(decimal.doubleValue())
      case unsupported => throw new UnsupportedOperationException(s"Cannot deserialize $unsupported to a WomValue")
    }
  }
}

object WomSingleFileJsonFormatter extends DefaultJsonProtocol {
  implicit object WomSingleFileJsonFormat extends RootJsonFormat[WomSingleFile] {
    def write(value: WomSingleFile) = JsString(value.value)
    def read(value: JsValue): WomSingleFile = value match {
      case JsString(path) => WomSingleFile(path)
      case unsupported => throw new UnsupportedOperationException(s"Cannot deserialize $unsupported to a WdlFile")
    }
  }
}

