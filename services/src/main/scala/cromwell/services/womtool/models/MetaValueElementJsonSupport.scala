package cromwell.services.womtool.models

import io.circe.{Encoder, Json}

import wom.callable.MetaValueElement
import wom.callable.MetaValueElement._

object MetaValueElementJsonSupport {

  // We only implement the encoder, because currently, limited query functionality is required
  // for the meta sections. Even in the encoder, we only support the type, and not the actual value.
  implicit val metaValueElementEncoder: Encoder[MetaValueElement] = new Encoder[MetaValueElement] {
    final def apply(a : MetaValueElement) : Json = {
      a match {
        case MetaValueElementNull =>
          Json.Null

        case MetaValueElementBoolean(b) =>
          Json.fromBoolean(b)

        case MetaValueElementFloat(x) =>
          Json.fromDouble(x).get

        case MetaValueElementInteger(a) =>
          Json.fromInt(a)

        case MetaValueElementString(s) =>
          Json.fromString(s)

        case MetaValueElementArray(vec) =>
          Json.fromValues(vec.map(apply(_)))

        case MetaValueElementObject(m) =>
          Json.fromFields(m.map{ case (key, value) =>
                            key -> apply(value)
                          })
      }
    }
  }
}
