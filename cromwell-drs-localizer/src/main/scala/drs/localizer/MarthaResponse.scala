package drs.localizer

import io.circe.{Decoder, Json, JsonObject}
import io.circe.generic.semiauto.deriveDecoder


object MarthaResponseJsonSupport {
  implicit val googleServiceAccountFormat: Decoder[GoogleServiceAccount] = deriveDecoder
  implicit val marthaResponseFormat: Decoder[MarthaResponse] = deriveDecoder
}

case class GoogleServiceAccount(data: Json)
case class MarthaResponse(contentType: Option[String],
                          size: Option[Long],
                          timeCreated: Option[String],
                          timeUpdated: Option[String],
                          bucket: Option[String],
                          name: Option[String],
                          gsUri: Option[String],
                          googleServiceAccount: Option[GoogleServiceAccount],
                          hashes: Option[JsonObject])

