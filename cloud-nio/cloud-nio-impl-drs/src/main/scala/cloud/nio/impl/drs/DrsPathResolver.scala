package cloud.nio.impl.drs

import cats.data.NonEmptyList
import cats.effect.{IO, Resource}
import cats.implicits._
import cloud.nio.impl.drs.MarthaResponseSupport._
import common.exception.toIO
import io.circe._
import io.circe.generic.semiauto._
import io.circe.parser.decode
import io.circe.syntax._
import mouse.boolean._
import org.apache.commons.lang3.exception.ExceptionUtils
import org.apache.http.client.methods.HttpPost
import org.apache.http.entity.{ContentType, StringEntity}
import org.apache.http.impl.client.HttpClientBuilder
import org.apache.http.util.EntityUtils
import org.apache.http.{HttpResponse, HttpStatus, StatusLine}

abstract class DrsPathResolver(drsConfig: DrsConfig, httpClientBuilder: HttpClientBuilder) {

  def getAccessToken: String

  private def makeHttpRequestToMartha(drsPath: String, fields: NonEmptyList[MarthaField.Value]): HttpPost = {
    val postRequest = new HttpPost(drsConfig.marthaUrl)
    val requestJson = MarthaRequest(drsPath, fields).asJson.noSpaces
    postRequest.setEntity(new StringEntity(requestJson, ContentType.APPLICATION_JSON))
    postRequest.setHeader("Authorization", s"Bearer $getAccessToken")
    postRequest
  }

  private def httpResponseToMarthaResponse(drsPathForDebugging: String)(httpResponse: HttpResponse): IO[MarthaResponse] = {
    val marthaResponseEntityOption = Option(httpResponse.getEntity).map(EntityUtils.toString)
    val responseStatusLine = httpResponse.getStatusLine

    val exceptionMsg = errorMessageFromResponse(drsPathForDebugging, marthaResponseEntityOption, responseStatusLine, drsConfig.marthaUrl)
    val responseEntityOption = (responseStatusLine.getStatusCode == HttpStatus.SC_OK).valueOrZero(marthaResponseEntityOption)
    val responseContentIO = toIO(responseEntityOption, exceptionMsg)

    responseContentIO.flatMap{ responseContent =>
      IO.fromEither(decode[MarthaResponse](responseContent))
    }.handleErrorWith {
      e => IO.raiseError(new RuntimeException(s"Unexpected response during DRS resolution: ${ExceptionUtils.getMessage(e)}"))
    }
  }

  private def executeMarthaRequest(httpPost: HttpPost): Resource[IO, HttpResponse]= {
    for {
      httpClient <- Resource.fromAutoCloseable(IO(httpClientBuilder.build()))
      httpResponse <- Resource.fromAutoCloseable(IO(httpClient.execute(httpPost)))
    } yield httpResponse
  }

  def rawMarthaResponse(drsPath: String, fields: NonEmptyList[MarthaField.Value]): Resource[IO, HttpResponse] = {
    val httpPost = makeHttpRequestToMartha(drsPath, fields)
    executeMarthaRequest(httpPost)
  }

  /** *
    * Resolves the DRS path through Martha url provided in the config.
    * Please note, this method returns an IO that would make a synchronous HTTP request to Martha when run.
    */
  def resolveDrsThroughMartha(drsPath: String, fields: NonEmptyList[MarthaField.Value]): IO[MarthaResponse] = {
    rawMarthaResponse(drsPath, fields).use(httpResponseToMarthaResponse(drsPathForDebugging = drsPath))
  }
}

final case class DrsConfig(marthaUrl: String)

object MarthaField extends Enumeration {
  val GsUri: MarthaField.Value = Value("gsUri")
  val Size: MarthaField.Value = Value("size")
  val TimeCreated: MarthaField.Value = Value("timeCreated")
  val TimeUpdated: MarthaField.Value = Value("timeUpdated")
  val GoogleServiceAccount: MarthaField.Value = Value("googleServiceAccount")
  val Hashes: MarthaField.Value = Value("hashes")
  val FileName: MarthaField.Value = Value("fileName")
}

final case class MarthaRequest(url: String, fields: NonEmptyList[MarthaField.Value])

final case class Url(url: String)
final case class ChecksumObject(checksum: String, `type`: String)
final case class DosDataObject(name: Option[String],
                               size: Option[Long],
                               checksums: Option[Array[ChecksumObject]],
                               created: Option[String],
                               updated: Option[String],
                               urls: Array[Url])
final case class DosObject(data_object: DosDataObject)
final case class SADataObject(data: Json)

final case class MarthaV2Response(dos: DosObject, googleServiceAccount: Option[SADataObject])

/**
  * A response from `martha_v3` or converted from `martha_v2`.
  *
  * @param size Size of the object stored at gsUri
  * @param timeCreated The creation time of the object at gsUri
  * @param timeUpdated The last update time of the object at gsUri
  * @param gsUri Where the object bytes are stored, possibly using a generated path name such as "gs://bucket/12/345"
  * @param googleServiceAccount The service account to access the gsUri contents
  * @param fileName A possible different file name for the object at gsUri, ex: "gsutil cp gs://bucket/12/345 my.vcf"
  * @param hashes Hashes for the contents stored at gsUri
  */
final case class MarthaResponse(size: Option[Long],
                                timeCreated: Option[String],
                                timeUpdated: Option[String],
                                gsUri: Option[String],
                                googleServiceAccount: Option[SADataObject],
                                fileName: Option[String],
                                hashes: Option[Map[String, String]])

// Adapted from https://github.com/broadinstitute/martha/blob/f31933a3a11e20d30698ec4b4dc1e0abbb31a8bc/common/helpers.js#L210-L218
final case class MarthaFailureResponse(response: MarthaFailureResponsePayload)
final case class MarthaFailureResponsePayload(text: String)

object MarthaResponseSupport {

  implicit lazy val marthaFieldEncoder: Encoder[MarthaField.Value] = Encoder.encodeEnumeration(MarthaField)
  implicit lazy val marthaRequestEncoder: Encoder[MarthaRequest] = deriveEncoder

  implicit lazy val urlDecoder: Decoder[Url] = deriveDecoder
  implicit lazy val checksumObjectDecoder: Decoder[ChecksumObject] = deriveDecoder
  implicit lazy val dosDataObjectDecoder: Decoder[DosDataObject] = deriveDecoder
  implicit lazy val dosObjectDecoder: Decoder[DosObject] = deriveDecoder
  implicit lazy val saDataObjectDecoder: Decoder[SADataObject] = deriveDecoder
  private lazy val marthaV3ResponseDecoder: Decoder[MarthaResponse] = deriveDecoder
  private lazy val marthaV2ResponseDecoder: Decoder[MarthaResponse] =
    deriveDecoder[MarthaV2Response] map convertMarthaResponseV2ToV3
  implicit lazy val marthaResponseDecoder: Decoder[MarthaResponse] =
    marthaV2ResponseDecoder or marthaV3ResponseDecoder

  implicit lazy val marthaFailureResponseDecoder: Decoder[MarthaFailureResponse] = deriveDecoder
  implicit lazy val marthaFailureResponsePayloadDecoder: Decoder[MarthaFailureResponsePayload] = deriveDecoder

  private val GcsScheme = "gs://"

  private def convertChecksumsToHashesMap(checksums: Array[ChecksumObject]): Map[String, String] = {
    checksums.flatMap (checksumObj => Map(checksumObj.`type` -> checksumObj.checksum)).toMap
  }

  def getGcsBucketAndName(gcsUrl: String): (String, String) = {
     val array = gcsUrl.substring(GcsScheme.length).split("/", 2)
      (array(0), array(1))
  }

  def convertMarthaResponseV2ToV3(response: MarthaV2Response): MarthaResponse = {
    val dataObject = response.dos.data_object
    val fileName = dataObject.name
    val size = dataObject.size
    val timeCreated = dataObject.created
    val timeUpdated = dataObject.updated
    val hashesMap = dataObject.checksums.map(convertChecksumsToHashesMap)
    val gcsUrl = dataObject.urls.find(_.url.startsWith(GcsScheme)).map(_.url)

    MarthaResponse(
      size = size,
      timeCreated = timeCreated,
      timeUpdated = timeUpdated,
      gsUri = gcsUrl,
      googleServiceAccount = response.googleServiceAccount,
      fileName = fileName,
      hashes = hashesMap
    )
  }

  def errorMessageFromResponse(drsPathForDebugging: String, marthaResponseEntityOption: Option[String], responseStatusLine: StatusLine, marthaUri: String): String = {
    val baseMessage = s"Could not access object \'$drsPathForDebugging\'. Status: ${responseStatusLine.getStatusCode}, reason: \'${responseStatusLine.getReasonPhrase}\', Martha location: \'$marthaUri\', message: "

    marthaResponseEntityOption match {
      case Some(entity) =>
        val maybeErrorResponse: Either[Error, MarthaFailureResponse] = decode[MarthaFailureResponse](entity)
        maybeErrorResponse match {
          case Left(_) =>
            // Not parsable as a `MarthaFailureResponse`
            baseMessage + s"\'$entity\'"
          case Right(decoded) =>
            // Is a `MarthaFailureResponse`
            baseMessage + s"\'${decoded.response.text}\'"
        }
      case None =>
        // No entity in HTTP response
        baseMessage + "(empty response)"
    }
  }
}
