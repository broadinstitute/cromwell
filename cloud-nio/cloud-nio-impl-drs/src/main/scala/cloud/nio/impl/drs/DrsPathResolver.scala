package cloud.nio.impl.drs

import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import cats.effect.{IO, Resource}
import cats.implicits._
import cloud.nio.impl.drs.DrsPathResolver.{FatalRetryDisposition, RegularRetryDisposition}
import cloud.nio.impl.drs.DrsResolverResponseSupport._
import common.exception.{toIO, AggregatedMessageException}
import common.validation.ErrorOr.ErrorOr
import io.circe._
import io.circe.generic.semiauto._
import io.circe.parser.decode
import io.circe.syntax._
import mouse.boolean._
import org.apache.commons.lang3.exception.ExceptionUtils
import org.apache.http.client.methods.{HttpGet, HttpPost}
import org.apache.http.entity.{ContentType, StringEntity}
import org.apache.http.impl.client.HttpClientBuilder
import org.apache.http.impl.conn.PoolingHttpClientConnectionManager
import org.apache.http.util.EntityUtils
import org.apache.http.{HttpResponse, HttpStatus, StatusLine}

import java.nio.ByteBuffer
import java.nio.channels.{Channels, ReadableByteChannel}
import scala.util.Try

class DrsPathResolver(drsConfig: DrsConfig, drsCredentials: DrsCredentials) {

  protected lazy val httpClientBuilder: HttpClientBuilder = {
    val clientBuilder = HttpClientBuilder.create()
    val retryHandler = new DrsResolverHttpRequestRetryStrategy(drsConfig)
    clientBuilder
      .setRetryHandler(retryHandler)
      .setServiceUnavailableRetryStrategy(retryHandler)
    clientBuilder.setConnectionManager(connectionManager)
    clientBuilder.setConnectionManagerShared(true)
    clientBuilder
  }

  def getAccessToken: ErrorOr[String] = drsCredentials.getAccessToken

  private lazy val currentCloudPlatform: Option[DrsCloudPlatform.Value] = drsCredentials match {
    case _: GoogleOauthDrsCredentials => Option(DrsCloudPlatform.GoogleStorage)
    case GoogleAppDefaultTokenStrategy => Option(DrsCloudPlatform.GoogleStorage)
    case _: AzureDrsCredentials => Option(DrsCloudPlatform.Azure)
    case _ => None
  }

  def makeDrsResolverRequest(drsPath: String, fields: NonEmptyList[DrsResolverField.Value]): DrsResolverRequest =
    DrsResolverRequest(drsPath, currentCloudPlatform, fields)

  private def makeHttpRequestToDrsResolver(drsPath: String,
                                           fields: NonEmptyList[DrsResolverField.Value]
  ): Resource[IO, HttpPost] = {
    val io = getAccessToken match {
      case Valid(token) =>
        IO {
          val postRequest = new HttpPost(drsConfig.drsResolverUrl)
          val requestJson = makeDrsResolverRequest(drsPath, fields).asJson.noSpaces
          postRequest.setEntity(new StringEntity(requestJson, ContentType.APPLICATION_JSON))
          postRequest.setHeader("Authorization", s"Bearer $token")
          postRequest
        }
      case Invalid(errors) =>
        IO.raiseError(AggregatedMessageException("Error getting access token", errors.toList))
    }
    Resource.eval(io)
  }

  private def httpResponseToDrsResolverResponse(
    drsPathForDebugging: String
  )(httpResponse: HttpResponse): IO[DrsResolverResponse] = {
    val responseStatusLine = httpResponse.getStatusLine
    val status = responseStatusLine.getStatusCode

    lazy val retryMessage = {
      val reason = responseStatusLine.getReasonPhrase
      s"Unexpected response during resolution of '$drsPathForDebugging': HTTP status $status: $reason"
    }

    status match {
      case 408 | 429 =>
        IO.raiseError(new RuntimeException(retryMessage) with RegularRetryDisposition)
      case s if s / 100 == 4 =>
        IO.raiseError(new RuntimeException(retryMessage) with FatalRetryDisposition)
      case s if s / 100 == 5 =>
        IO.raiseError(new RuntimeException(retryMessage) with RegularRetryDisposition)
      case _ =>
        val drsResolverResponseEntityOption = Option(httpResponse.getEntity).map(EntityUtils.toString)
        val exceptionMsg = errorMessageFromResponse(drsPathForDebugging,
                                                    drsResolverResponseEntityOption,
                                                    responseStatusLine,
                                                    drsConfig.drsResolverUrl
        )
        val responseEntityOption =
          (responseStatusLine.getStatusCode == HttpStatus.SC_OK).valueOrZero(drsResolverResponseEntityOption)
        val responseContentIO = toIO(responseEntityOption, exceptionMsg)

        responseContentIO
          .flatMap { responseContent =>
            IO.fromEither(decode[DrsResolverResponse](responseContent))
          }
          .handleErrorWith { e =>
            IO.raiseError(
              new RuntimeException(s"Unexpected response during DRS resolution: ${ExceptionUtils.getMessage(e)}")
            )
          }
    }
  }

  private def executeDrsResolverRequest(httpPost: HttpPost): Resource[IO, HttpResponse] =
    for {
      httpClient <- Resource.fromAutoCloseable(IO(httpClientBuilder.build()))
      httpResponse <- Resource.fromAutoCloseable(IO(httpClient.execute(httpPost)))
    } yield httpResponse

  def rawDrsResolverResponse(drsPath: String,
                             fields: NonEmptyList[DrsResolverField.Value]
  ): Resource[IO, HttpResponse] =
    for {
      httpPost <- makeHttpRequestToDrsResolver(drsPath, fields)
      response <- executeDrsResolverRequest(httpPost)
    } yield response

  /** *
    * Resolves the DRS path through DRS Resolver url provided in the config.
    * Please note, this method returns an IO that would make a synchronous HTTP request to DRS Resolver when run.
    */
  def resolveDrs(drsPath: String, fields: NonEmptyList[DrsResolverField.Value]): IO[DrsResolverResponse] =
    rawDrsResolverResponse(drsPath, fields).use(
      httpResponseToDrsResolverResponse(drsPathForDebugging = drsPath)
    )

  def openChannel(accessUrl: AccessUrl): IO[ReadableByteChannel] =
    IO {
      val httpGet = new HttpGet(accessUrl.url)
      accessUrl.headers.getOrElse(Map.empty).toList foreach { case (name, value) =>
        httpGet.addHeader(name, value)
      }
      val client = httpClientBuilder.build()
      val response = client.execute(httpGet)
      val inner = Channels.newChannel(response.getEntity.getContent)
      /*
      Create a wrapper ReadableByteChannel. When .close() is invoked on the wrapper, the wrapper will internally call
      .close() on:
      - the inner ReadableByteChannel
      - the HttpResponse
      - the HttpClient

      This ensures that when the channel is released any underlying HTTP connections are also released.
       */
      new ReadableByteChannel {
        override def read(dst: ByteBuffer): Int = inner.read(dst)

        override def isOpen: Boolean = inner.isOpen

        // noinspection ScalaUnusedExpression
        override def close(): Unit = {
          val innerTry = Try(inner.close())
          val responseTry = Try(response.close())
          val clientTry = Try(client.close())
          innerTry.get
          responseTry.get
          clientTry.get
        }
      }
    }
}

object DrsPathResolver {
  final val ExtractUriErrorMsg = "No access URL nor GCS URI starting with 'gs://' found in the DRS Resolver response!"
  sealed trait RetryDisposition
  // Should immediately fail the download attempt.
  trait FatalRetryDisposition extends RetryDisposition
  // Should increase the attempt counter and continue retrying if more retry attempts remain.
  trait RegularRetryDisposition extends RetryDisposition
}

object DrsResolverField extends Enumeration {
  val GsUri: DrsResolverField.Value = Value("gsUri")
  val Size: DrsResolverField.Value = Value("size")
  val TimeCreated: DrsResolverField.Value = Value("timeCreated")
  val TimeUpdated: DrsResolverField.Value = Value("timeUpdated")
  val BondProvider: DrsResolverField.Value = Value("bondProvider")
  val GoogleServiceAccount: DrsResolverField.Value = Value("googleServiceAccount")
  val Hashes: DrsResolverField.Value = Value("hashes")
  val FileName: DrsResolverField.Value = Value("fileName")
  val AccessUrl: DrsResolverField.Value = Value("accessUrl")
  val LocalizationPath: DrsResolverField.Value = Value("localizationPath")
}

// We supply a cloud platform value to the DRS service. In cases where the DRS repository
// has multiple cloud files associated with a DRS link, it will prefer sending a file on the same
// platform as this Cromwell instance. That is, if a DRS file has copies on both GCP and Azure,
// we'll get the GCP one when running on GCP and the Azure one when running on Azure.
object DrsCloudPlatform extends Enumeration {
  val GoogleStorage: DrsCloudPlatform.Value = Value("gs")
  val Azure: DrsCloudPlatform.Value = Value("azure")
  val AmazonS3: DrsCloudPlatform.Value = Value("s3") // supported by DRSHub but not currently used by us
}

final case class DrsResolverRequest(url: String,
                                    cloudPlatform: Option[DrsCloudPlatform.Value],
                                    fields: NonEmptyList[DrsResolverField.Value]
)

final case class SADataObject(data: Json)

final case class AccessUrl(url: String, headers: Option[Map[String, String]])

/**
  * A response from `drshub_v4`.
  *
  * @param size Size of the object stored at gsUri
  * @param timeCreated The creation time of the object at gsUri
  * @param timeUpdated The last update time of the object at gsUri
  * @param gsUri Where the object bytes are stored, possibly using a generated path name such as "gs://bucket/12/345"
  * @param bondProvider The bond provider returning the googleServiceAccount
  * @param googleServiceAccount The service account to access the gsUri contents created via bondProvider
  * @param fileName A possible different file name for the object at gsUri, ex: "gsutil cp gs://bucket/12/345 my.vcf"
  * @param hashes Hashes for the contents stored at gsUri
  * @param accessUrl URL to query for signed URL
  * @param localizationPath Optional localization path. TDR is currently the sole DRS provider specifying this value in
  *                         DRS metadata, via the `aliases` field. As this is a distinct field from `fileName` in DRS
  *                         metadata it is also made a distinct field in this response object.
  */
final case class DrsResolverResponse(size: Option[Long] = None,
                                     timeCreated: Option[String] = None,
                                     timeUpdated: Option[String] = None,
                                     gsUri: Option[String] = None,
                                     bondProvider: Option[String] = None,
                                     googleServiceAccount: Option[SADataObject] = None,
                                     fileName: Option[String] = None,
                                     hashes: Option[Map[String, String]] = None,
                                     accessUrl: Option[AccessUrl] = None,
                                     localizationPath: Option[String] = None
)

final case class DrsResolverFailureResponse(response: DrsResolverFailureResponsePayload)
final case class DrsResolverFailureResponsePayload(text: String)

object DrsResolverResponseSupport {

  implicit lazy val drsResolverFieldEncoder: Encoder[DrsResolverField.Value] =
    Encoder.encodeEnumeration(DrsResolverField)
  implicit lazy val drsResolverCloudPlatformEncoder: Encoder[DrsCloudPlatform.Value] =
    Encoder.encodeEnumeration(DrsCloudPlatform)
  implicit lazy val drsResolverRequestEncoder: Encoder[DrsResolverRequest] = deriveEncoder

  implicit lazy val saDataObjectDecoder: Decoder[SADataObject] = deriveDecoder
  implicit lazy val drsResolverResponseDecoder: Decoder[DrsResolverResponse] = deriveDecoder

  implicit lazy val drsResolverFailureResponseDecoder: Decoder[DrsResolverFailureResponse] = deriveDecoder
  implicit lazy val drsResolverFailureResponsePayloadDecoder: Decoder[DrsResolverFailureResponsePayload] = deriveDecoder

  implicit lazy val drsResolverAccessUrlDecoder: Decoder[AccessUrl] = deriveDecoder

  private val GcsScheme = "gs://"

  def getGcsBucketAndName(gcsUrl: String): (String, String) = {
    val array = gcsUrl.substring(GcsScheme.length).split("/", 2)
    (array(0), array(1))
  }

  def errorMessageFromResponse(drsPathForDebugging: String,
                               drsResolverResponseEntityOption: Option[String],
                               responseStatusLine: StatusLine,
                               drsResolverUri: String
  ): String = {
    val baseMessage =
      s"Could not access object \'$drsPathForDebugging\'. Status: ${responseStatusLine.getStatusCode}, reason: \'${responseStatusLine.getReasonPhrase}\', DRS Resolver location: \'$drsResolverUri\', message: "

    drsResolverResponseEntityOption match {
      case Some(entity) =>
        val maybeErrorResponse: Either[Error, DrsResolverFailureResponse] = decode[DrsResolverFailureResponse](entity)
        maybeErrorResponse match {
          case Left(_) =>
            // Not parsable as a `DrsResolverFailureResponse`
            baseMessage + s"\'$entity\'"
          case Right(decoded) =>
            // Is a `DrsResolverFailureResponse`
            baseMessage + s"\'${decoded.response.text}\'"
        }
      case None =>
        // No entity in HTTP response
        baseMessage + "(empty response)"
    }
  }

  lazy val connectionManager = {
    val connManager = new PoolingHttpClientConnectionManager()
    connManager.setMaxTotal(250)
    // Since the HttpClient is always talking to DRSHub,
    // make the max connections per route the same as max total connections
    connManager.setDefaultMaxPerRoute(250)
    connManager
  }
}
