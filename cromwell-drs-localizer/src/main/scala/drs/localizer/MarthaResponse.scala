package drs.localizer

import io.circe.{Decoder, Json}
import io.circe.generic.semiauto.deriveDecoder


object MarthaResponseJsonSupport {
  implicit val urlFormat: Decoder[Url] = deriveDecoder
  implicit val dataObject: Decoder[DrsDataObject] = deriveDecoder
  implicit val drsObjectFormat: Decoder[DrsObject] = deriveDecoder
  implicit val googleServiceAccountFormat: Decoder[GoogleServiceAccount] = deriveDecoder
  implicit val marthaResponseFormat: Decoder[MarthaResponse] = Decoder.forProduct2("dos", "googleServiceAccount")(MarthaResponse.apply)

  implicit val samErrorResponseFormat: Decoder[SamErrorResponse] = deriveDecoder
  implicit val samErrorResponseCodeFormat: Decoder[SamErrorResponseCode] = deriveDecoder
  implicit val marthaErrorResponseFormat: Decoder[MarthaErrorResponse] = deriveDecoder
}

case class Url(url: String)
case class DrsDataObject(urls: Array[Url])
case class DrsObject(data_object: DrsDataObject)
case class GoogleServiceAccount(data: Json)
case class MarthaResponse(drs: DrsObject, googleServiceAccount: Option[GoogleServiceAccount])


case class SamErrorResponse(text: String)
case class SamErrorResponseCode(status: Int)
case class MarthaErrorResponse(response: SamErrorResponse)

