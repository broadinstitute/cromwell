package cwl

import spray.json.{JsArray, JsBoolean, JsNull, JsNumber, JsObject, JsString, JsValue, RootJsonFormat}
import wom.values.{WomArray, WomBoolean, WomFile, WomFloat, WomInteger, WomMap, WomObjectLike, WomOptionalValue, WomPair, WomString, WomValue}

// Sometimes it's hard for compiler plugins to find this in ScatterLogicSpec, so here it is again:
object WomValueJsonFormatter {
  implicit object WomValueJsonFormat extends RootJsonFormat[WomValue] {
    def write(value: WomValue): JsValue = value match {
      case s: WomString => JsString(s.value)
      case i: WomInteger => JsNumber(i.value)
      case f: WomFloat => JsNumber(f.value)
      case b: WomBoolean => JsBoolean(b.value)
      case f: WomFile => JsString(f.value)
      case o: WomObjectLike => new JsObject(o.values map {case(k, v) => k -> write(v)})
      case a: WomArray => new JsArray(a.value.map(write).toVector)
      case m: WomMap => new JsObject(m.value map {case(k,v) => k.valueString -> write(v)})
      case q: WomPair => new JsObject(Map("left" -> write(q.left), "right" -> write(q.right)))
      case WomOptionalValue(_, Some(innerValue)) => write(innerValue)
      case WomOptionalValue(_, None) => JsNull
      // handles WdlExpression
      case v: WomValue => JsString(v.toWomString)

    }
    def read(value: JsValue): WomValue = throw new UnsupportedOperationException
  }
}
