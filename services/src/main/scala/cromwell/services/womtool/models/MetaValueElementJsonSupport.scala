package cromwell.services.womtool.models

import spray.json._
import wom.callable.MetaValueElement
import wom.callable.MetaValueElement._

object MetaValueElementJsonSupport {

  // We only implement the writer, because currently, the reader is not required.
  implicit val metaValueElementEncoder: JsonWriter[MetaValueElement] = new JsonWriter[MetaValueElement] {
    final def write(a: MetaValueElement): JsValue =
      a match {
        case MetaValueElementNull =>
          JsNull

        case MetaValueElementBoolean(b) =>
          JsBoolean(b)

        case MetaValueElementFloat(x) =>
          JsNumber(x)

        case MetaValueElementInteger(a) =>
          JsNumber(a)

        case MetaValueElementString(s) =>
          JsString(s)

        case MetaValueElementArray(vec) =>
          JsArray(vec.map(write))

        case MetaValueElementObject(m) =>
          JsObject(m.map { case (key, value) =>
            key -> write(value)
          })
      }
  }
}
