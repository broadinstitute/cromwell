package drs.localizer

import io.circe.{Decoder, Json}
import io.circe.generic.semiauto.deriveDecoder


object MarthaResponseJsonSupport {
  implicit val urlFormat: Decoder[Url] = deriveDecoder
  implicit val dataObject: Decoder[DrsDataObject] = deriveDecoder
  implicit val drsObjectFormat: Decoder[DrsObject] = deriveDecoder
  implicit val googleServiceAccountFormat: Decoder[GoogleServiceAccount] = deriveDecoder
  // Martha is still returning objects keyed by the obsolete "dos" terminology rather than the current term "drs".
  // In order to avoid having Cromwell's case classes use the obsolete terminology that would arise from a derived
  // decoder, this `forProduct2` construct instructs Circe to take the value keyed by `dos` and pass that as the
  // first argument to `MarthaResponse.apply`, which happens to be the constructor parameter formally named `drs`.
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

