package drs.localizer

import io.circe.{Decoder, Json}
import io.circe.generic.semiauto.deriveDecoder


object MarthaResponseJsonSupport {
  implicit val urlFormat: Decoder[Url] = deriveDecoder
  implicit val dataObject: Decoder[DosDataObject] = deriveDecoder
  implicit val dosObjectFormat: Decoder[DosObject] = deriveDecoder
  implicit val googleServiceAccountFormat: Decoder[GoogleServiceAccount] = deriveDecoder
  implicit val marthaResponseFormat: Decoder[MarthaResponse] = deriveDecoder

  implicit val samErrorResponseFormat: Decoder[SamErrorResponse] = deriveDecoder
  implicit val samErrorResponseCodeFormat: Decoder[SamErrorResponseCode] = deriveDecoder
  implicit val marthaErrorResponseFormat: Decoder[MarthaErrorResponse] = deriveDecoder
}

case class Url(url: String)
case class DosDataObject(urls: Array[Url])
case class DosObject(data_object: DosDataObject)
case class GoogleServiceAccount(data: Json)
case class MarthaResponse(dos: DosObject, googleServiceAccount: Option[GoogleServiceAccount])


case class SamErrorResponse(text: String)
case class SamErrorResponseCode(status: Int)
case class MarthaErrorResponse(response: SamErrorResponse)

