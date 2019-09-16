package cloud.nio.impl.drs

import cats.effect.{IO, Resource}
import cats.instances.option._
import cats.instances.string._
import common.exception._
import io.circe._
import io.circe.generic.semiauto._
import io.circe.parser.decode
import mouse.boolean._
import org.apache.commons.lang3.exception.ExceptionUtils
import org.apache.http.client.methods.HttpPost
import org.apache.http.entity.{ContentType, StringEntity}
import org.apache.http.impl.client.HttpClientBuilder
import org.apache.http.util.EntityUtils
import org.apache.http.{HttpResponse, HttpStatus}

case class DrsPathResolver(drsConfig: DrsConfig, httpClientBuilder: HttpClientBuilder) {

  implicit lazy val urlDecoder: Decoder[Url] = deriveDecoder
  implicit lazy val checksumDecoder: Decoder[ChecksumObject] = deriveDecoder
  implicit lazy val dataObjectDecoder: Decoder[DrsDataObject] = deriveDecoder
  implicit lazy val drsObjectDecoder: Decoder[DrsObject] = deriveDecoder
  implicit lazy val saDataObjectDecoder: Decoder[SADataObject] = deriveDecoder
  // Martha is still returning objects keyed by the obsolete "dos" terminology rather than the current term "drs".
  // In order to avoid having Cromwell's case classes use the obsolete terminology that would arise from a derived
  // decoder, this `forProduct2` construct instructs Circe to take the value keyed by `dos` and pass that as the
  // first argument to `MarthaResponse.apply`, which happens to be the constructor parameter formally named `drs`.
  implicit lazy val marthaResponseDecoder: Decoder[MarthaResponse] = Decoder.forProduct2("dos", "googleServiceAccount")(MarthaResponse.apply)

  private val DrsPathToken = "${drsPath}"


  private def makeHttpRequestToMartha(drsPath: String, serviceAccount: Option[String]): HttpPost = {
    val postRequest = new HttpPost(drsConfig.marthaUri)
    val requestJson = drsConfig.marthaRequestJsonTemplate.replace(DrsPathToken, drsPath)
    postRequest.setEntity(new StringEntity(requestJson, ContentType.APPLICATION_JSON))
    serviceAccount.foreach(sa => postRequest.setHeader("Authorization", s"Bearer $sa"))
    postRequest
  }


  private def httpResponseToMarthaResponse(httpResponse: HttpResponse): IO[MarthaResponse] = {
    val marthaResponseEntityOption = Option(httpResponse.getEntity).map(EntityUtils.toString)
    val responseStatusLine = httpResponse.getStatusLine

    val exceptionMsg = s"Unexpected response resolving DRS path through Martha url ${drsConfig.marthaUri}. Error: ${responseStatusLine.getStatusCode} ${responseStatusLine.getReasonPhrase}."
    val responseEntityOption = (responseStatusLine.getStatusCode == HttpStatus.SC_OK).valueOrZero(marthaResponseEntityOption)
    val responseContentIO = toIO(responseEntityOption, exceptionMsg)

    responseContentIO.flatMap(responseContent => IO.fromEither(decode[MarthaResponse](responseContent))).handleErrorWith {
      e => IO.raiseError(new RuntimeException(s"Failed to parse response from Martha into a case class. Error: ${ExceptionUtils.getMessage(e)}"))
    }
  }

  private def executeMarthaRequest(httpPost: HttpPost): Resource[IO, HttpResponse]= {
    for {
      httpClient <- Resource.fromAutoCloseable(IO(httpClientBuilder.build()))
      httpResponse <- Resource.fromAutoCloseable(IO(httpClient.execute(httpPost)))
    } yield httpResponse
  }


  def rawMarthaResponse(drsPath: String, serviceAccount: Option[String] = None): Resource[IO, HttpResponse] = {
    val httpPost = makeHttpRequestToMartha(drsPath, serviceAccount)
    executeMarthaRequest(httpPost)
  }

  /** *
    * Resolves the DRS path through Martha url provided in the config.
    * Please note, this method returns an IO that would make a synchronous HTTP request to Martha when run.
    */
  def resolveDrsThroughMartha(drsPath: String, serviceAccount: Option[String] = None): IO[MarthaResponse] = {
    rawMarthaResponse(drsPath, serviceAccount).use(httpResponseToMarthaResponse)
  }
}

case class DrsConfig(marthaUri: String, marthaRequestJsonTemplate: String)


case class Url(url: String)

case class ChecksumObject(checksum: String, `type`: String)

case class DrsDataObject(size: Option[Long],
                         checksums: Option[Array[ChecksumObject]],
                         updated: Option[String],
                         urls: Array[Url])

case class DrsObject(data_object: DrsDataObject)

case class SADataObject(data: Json)

case class MarthaResponse(drs: DrsObject, googleServiceAccount: Option[SADataObject])
