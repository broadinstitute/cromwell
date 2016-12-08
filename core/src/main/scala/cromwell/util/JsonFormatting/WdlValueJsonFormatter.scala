package cromwell.util.JsonFormatting

import spray.json._
import wdl4s.WdlExpression
import wdl4s.types.{WdlArrayType, WdlMapType, WdlStringType}
import wdl4s.values._

object WdlValueJsonFormatter extends DefaultJsonProtocol {
  implicit object WdlValueJsonFormat extends RootJsonFormat[WdlValue] {
    def write(value: WdlValue): JsValue = value match {
      case s: WdlString => JsString(s.value)
      case i: WdlInteger => JsNumber(i.value)
      case f: WdlFloat => JsNumber(f.value)
      case b: WdlBoolean => JsBoolean(b.value)
      case f: WdlFile => JsString(f.value)
      case o: WdlObject => new JsObject(o.value map {case(k, v) => k -> write(v)})
      case a: WdlArray => new JsArray(a.value.map(write).toVector)
      case m: WdlMap => new JsObject(m.value map {case(k,v) => k.valueString -> write(v)})
      case e: WdlExpression => JsString(e.toWdlString)
      case q: WdlPair => new JsObject(Map("left" -> write(q.left), "right" -> write(q.right)))
      case WdlOptionalValue(_, Some(innerValue)) => write(innerValue)
      case WdlOptionalValue(_, None) => JsNull
    }

    // NOTE: This assumes a map's keys are strings. Since we're coming from JSON this is fine.
    // This won't support a map with complex keys (e.g. WdlMapType(WdlMapType(WdlIntegerType, WdlIntegerType), WdlStringType)
    // That would require a more inventive JSON which splits out the key and value as full fields in their own right...
    // In addition, we make a lot of assumptions about what type of WdlValue to create. Oh well... it should all fall out in the coercion (fingercrossed)!
    def read(value: JsValue): WdlValue = value match {
      case JsObject(fields) =>
        val wdlFields: Map[WdlValue, WdlValue] = fields map {case (k, v) => WdlString(k) -> read(v)}
        if (fields.isEmpty) WdlMap(WdlMapType(WdlStringType, WdlStringType), Map.empty[WdlValue, WdlValue])
        else WdlMap(WdlMapType(wdlFields.head._1.wdlType, wdlFields.head._2.wdlType), wdlFields)
      case JsArray(vector) if vector.nonEmpty => WdlArray(WdlArrayType(read(vector.head).wdlType), vector map read)
      case JsString(str) => WdlString(str)
      case JsBoolean(bool) => WdlBoolean(bool)
      case JsNumber(decimal) if decimal.isValidInt => WdlInteger(decimal.toIntExact)
      case JsNumber(decimal) => WdlFloat(decimal.doubleValue())
      case unsupported => throw new UnsupportedOperationException(s"Cannot deserialize $unsupported to a WdlValue")
    }
  }
}

object WdlFileJsonFormatter extends DefaultJsonProtocol {
  implicit object WdlFileJsonFormat extends RootJsonFormat[WdlFile] {
    def write(value: WdlFile) = JsString(value.value)
    def read(value: JsValue): WdlFile = value match {
      case JsString(path) => WdlFile(path)
      case unsupported => throw new UnsupportedOperationException(s"Cannot deserialize $unsupported to a WdlFile")
    }
  }
}

