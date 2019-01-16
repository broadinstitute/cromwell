package cromwell.services.womtool.models

import io.circe.{Decoder, Encoder, HCursor, Json}
import wom.types.{WomStringType, WomType}

object WomTypeJsonSupport {
  implicit val womTypeEncoder: Encoder[WomType] = new Encoder[WomType] {
    final def apply(a: WomType): Json = Json.obj(
      ("typeName", Json.fromString(a.toDisplayString))
    )
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
