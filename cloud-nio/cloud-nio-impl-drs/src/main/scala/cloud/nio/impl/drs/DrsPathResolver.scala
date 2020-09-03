package cloud.nio.impl.drs

import cats.effect.{IO, Resource}
import cloud.nio.impl.drs.MarthaResponseSupport._
import common.exception.toIO
import io.circe._
import io.circe.generic.semiauto._
import io.circe.parser.decode
import mouse.boolean._
import org.apache.commons.lang3.exception.ExceptionUtils
import org.apache.http.client.methods.HttpPost
import org.apache.http.entity.{ContentType, StringEntity}
import org.apache.http.impl.client.HttpClientBuilder
import org.apache.http.util.EntityUtils
import org.apache.http.{HttpResponse, HttpStatus, StatusLine}


abstract class DrsPathResolver(drsConfig: DrsConfig, httpClientBuilder: HttpClientBuilder) {

  private val DrsPathToken = "${drsPath}"

  def getAccessToken: String

  private def makeHttpRequestToMartha(drsPath: String): HttpPost = {
    val postRequest = new HttpPost(drsConfig.marthaUri)
    val requestJson = drsConfig.marthaRequestJsonTemplate.replace(DrsPathToken, drsPath)
    postRequest.setEntity(new StringEntity(requestJson, ContentType.APPLICATION_JSON))
    postRequest.setHeader("Authorization", s"Bearer $getAccessToken")
    postRequest
  }

  private def httpResponseToMarthaResponse(drsPathForDebugging: String)(httpResponse: HttpResponse): IO[MarthaResponse] = {
    val marthaResponseEntityOption = Option(httpResponse.getEntity).map(EntityUtils.toString)
    val responseStatusLine = httpResponse.getStatusLine

    val exceptionMsg = errorMessageFromResponse(drsPathForDebugging, marthaResponseEntityOption, responseStatusLine, drsConfig.marthaUri)
    val responseEntityOption = (responseStatusLine.getStatusCode == HttpStatus.SC_OK).valueOrZero(marthaResponseEntityOption)
    val responseContentIO = toIO(responseEntityOption, exceptionMsg)

    responseContentIO.flatMap{ responseContent =>
      if (drsConfig.marthaUri.endsWith("martha_v3")) IO.fromEither(decode[MarthaResponse](responseContent))
      else IO.fromEither(decode[MarthaV2Response](responseContent).map(convertMarthaResponseV2ToV3))
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

  def rawMarthaResponse(drsPath: String): Resource[IO, HttpResponse] = {
    val httpPost = makeHttpRequestToMartha(drsPath)
    executeMarthaRequest(httpPost)
  }

  /** *
    * Resolves the DRS path through Martha url provided in the config.
    * Please note, this method returns an IO that would make a synchronous HTTP request to Martha when run.
    */
  def resolveDrsThroughMartha(drsPath: String): IO[MarthaResponse] = {
    rawMarthaResponse(drsPath).use(httpResponseToMarthaResponse(drsPathForDebugging = drsPath))
  }
}


final case class DrsConfig(marthaUri: String, marthaRequestJsonTemplate: String)

final case class Url(url: String)
final case class ChecksumObject(checksum: String, `type`: String)
final case class DrsDataObject(size: Option[Long],
                               checksums: Option[Array[ChecksumObject]],
                               updated: Option[String],
                               urls: Array[Url])
final case class DrsObject(data_object: DrsDataObject)
final case class SADataObject(data: Json)

final case class MarthaV2Response(dos: DrsObject, googleServiceAccount: Option[SADataObject])

final case class MarthaResponse(size: Option[Long],
                                timeUpdated: Option[String],
                                bucket: Option[String],
                                name: Option[String],
                                gsUri: Option[String],
                                googleServiceAccount: Option[SADataObject],
                                hashes: Option[Map[String, String]])

// Adapted from https://github.com/broadinstitute/martha/blob/f31933a3a11e20d30698ec4b4dc1e0abbb31a8bc/common/helpers.js#L210-L218
final case class MarthaFailureResponse(response: MarthaFailureResponsePayload)
final case class MarthaFailureResponsePayload(text: String)

object MarthaResponseSupport {

  implicit lazy val urlDecoder: Decoder[Url] = deriveDecoder
  implicit lazy val checksumDecoder: Decoder[ChecksumObject] = deriveDecoder
  implicit lazy val dataObjectDecoder: Decoder[DrsDataObject] = deriveDecoder
  implicit lazy val drsObjectDecoder: Decoder[DrsObject] = deriveDecoder
  implicit lazy val saDataObjectDecoder: Decoder[SADataObject] = deriveDecoder
  implicit lazy val marthaV2ResponseDecoder: Decoder[MarthaV2Response] = deriveDecoder
  implicit lazy val marthaResponseDecoder: Decoder[MarthaResponse] = deriveDecoder

  implicit lazy val marthaFailureResponseDecoder: Decoder[MarthaFailureResponse] = deriveDecoder
  implicit lazy val marthaFailureResponsePayloadDecoder: Decoder[MarthaFailureResponsePayload] = deriveDecoder

  private val GcsScheme = "gs://"

  private def convertChecksumsToHashesMap(checksums: Array[ChecksumObject]): Map[String, String] = {
    checksums.flatMap (checksumObj => Map(checksumObj.`type` -> checksumObj.checksum)).toMap
  }

  private def getGcsBucketAndName(gcsUrlOption: Option[String]): (Option[String], Option[String]) = {
    gcsUrlOption match {
      case Some(gcsUrl) =>
       val array = gcsUrl.substring(GcsScheme.length).split("/", 2)
        (Option(array(0)), Option(array(1)))
      case None => (None, None)
    }
  }

  def convertMarthaResponseV2ToV3(response: MarthaV2Response): MarthaResponse = {
    val dataObject = response.dos.data_object
    val size = dataObject.size
    val timeUpdated = dataObject.updated
    val hashesMap = dataObject.checksums.map(convertChecksumsToHashesMap)
    val gcsUrl = dataObject.urls.find(_.url.startsWith(GcsScheme)).map(_.url)
    val (bucketName, fileName) = getGcsBucketAndName(gcsUrl)

    MarthaResponse(
      size = size,
      timeUpdated = timeUpdated,
      bucket = bucketName,
      name = fileName,
      gsUri = gcsUrl,
      googleServiceAccount = response.googleServiceAccount,
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
