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
          Json.obj(
            ("typeName", Json.fromString("Null"))
          )

        case MetaValueElementBoolean(_) =>
          Json.obj(
            ("typeName", Json.fromString("Boolean"))
          )

        case MetaValueElementFloat(_) =>
          Json.obj(
            ("typeName", Json.fromString("Float"))
          )

        case MetaValueElementInteger(_) =>
          Json.obj(
            ("typeName", Json.fromString("Int"))
          )

        case MetaValueElementString(_) =>
          Json.obj(
            ("typeName", Json.fromString("String"))
          )

        case MetaValueElementArray(_) =>
          Json.obj(
            ("typeName", Json.fromString("Array"))
          )

        case MetaValueElementObject(_) =>
          Json.obj(
            ("typeName", Json.fromString("Object"))
          )
      }
    }
  }
}
