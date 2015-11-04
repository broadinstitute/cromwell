package cromwell.binding.values

import cromwell.binding.WdlExpression
import cromwell.binding.types.WdlArrayType
import spray.json._

object WdlValueJsonFormatter extends DefaultJsonProtocol {
  implicit object WdlValueJsonFormat extends RootJsonFormat[WdlValue] {
    def write(value: WdlValue) = value match {
      case s: WdlString => JsString(s.value)
      case i: WdlInteger => JsNumber(i.value)
      case f: WdlFloat => JsNumber(f.value)
      case b: WdlBoolean => JsBoolean(b.value)
      case f: WdlFile => JsString(f.value)
      case o: WdlObject => new JsObject(o.value map {case(k, v) => k -> write(v)})
      case a: WdlArray => new JsArray(a.value.map(write).toVector)
      case m: WdlMap => new JsObject(m.value map {case(k,v) => k.valueString -> write(v)})
      case e: WdlExpression => JsString(e.toWdlString)
    }
    // NOTE: This is NOT a completely safe way to read values from JSON. Should only be used for testing.
    def read(value: JsValue): WdlValue = value match {
      case JsObject(fields) => WdlObject(fields map {case (k, v) => k -> read(value)})
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
