package cromwell.services.womtool.models

import io.circe.{Decoder, Encoder, HCursor, Json}
import wom.types.{WomArrayType, WomOptionalType, WomStringType, WomType}

object WomTypeJsonSupport {
  implicit val womTypeEncoder: Encoder[WomType] = new Encoder[WomType] {
    final def apply(a: WomType): Json = {
      a match {
        case a: WomArrayType =>
          Json.obj(
            ("typeName", Json.fromString("Array")),
            ("arrayType", womTypeEncoder.apply(a.memberType))
          )
        case a: WomOptionalType =>
          Json.obj(
            ("typeName", Json.fromString("Optional")),
            ("optionalType", womTypeEncoder.apply(a.memberType))
          )
        case _ =>
          Json.obj(
            ("typeName", Json.fromString(a.toDisplayString))
          )
      }
    }
  }

  implicit val womTypeDecoder: Decoder[WomType] = new Decoder[WomType] {
    final def apply(c: HCursor): Decoder.Result[WomType] =
    for {
      _ <- c.downField("foo").as[String]
    } yield {
      WomStringType
    }
  }
}
