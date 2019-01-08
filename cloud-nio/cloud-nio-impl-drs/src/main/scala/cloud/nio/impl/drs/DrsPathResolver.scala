package cloud.nio.impl.drs

import cats.effect.{IO, Resource}
import io.circe._
import io.circe.generic.semiauto._
import io.circe.parser.decode
import org.apache.commons.lang3.exception.ExceptionUtils
import org.apache.http.client.methods.HttpPost
import org.apache.http.entity.{ContentType, StringEntity}
import org.apache.http.impl.client.HttpClientBuilder
import org.apache.http.util.EntityUtils
import org.apache.http.{HttpResponse, HttpStatus, StatusLine}

//Do not remove this import. If removed IntelliJ throws error for method 'executeMarthaRequest' stating
//'value map is not a member of cats.effect.Resource[cats.effect.IO,org.apache.http.client.methods.CloseableHttpResponse]'
import cats.syntax.functor._


case class DrsPathResolver(drsConfig: DrsConfig, httpClientBuilder: HttpClientBuilder) {

  implicit lazy val urlDecoder: Decoder[Url] = deriveDecoder
  implicit lazy val checksumDecoder: Decoder[ChecksumObject] = deriveDecoder
  implicit lazy val dataObjectDecoder: Decoder[DosDataObject] = deriveDecoder
  implicit lazy val dosObjectDecoder: Decoder[DosObject] = deriveDecoder
  implicit lazy val saDataObjectDecoder: Decoder[SADataObject] = deriveDecoder
  implicit lazy val marthaResponseDecoder: Decoder[MarthaResponse] = deriveDecoder

  private val DrsPathToken = "${drsPath}"

  private def unexpectedExceptionResponse(status: StatusLine): RuntimeException = {
    new RuntimeException(s"Unexpected response resolving DRS path through Martha url ${drsConfig.marthaUri}. Error: ${status.getStatusCode} ${status.getReasonPhrase}.")
  }


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

    val responseContentIO =
      IO.fromEither(
        marthaResponseEntityOption.filter { _ => responseStatusLine.getStatusCode == HttpStatus.SC_OK}
          .toRight(unexpectedExceptionResponse(responseStatusLine))
      )

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
    * Please note, this method makes a synchronous HTTP request to Martha.
    */
  def resolveDrsThroughMartha(drsPath: String, serviceAccount: Option[String] = None): IO[MarthaResponse] = {
    rawMarthaResponse(drsPath, serviceAccount).use(httpResponseToMarthaResponse)
  }
}

case class DrsConfig(marthaUri: String, marthaRequestJsonTemplate: String)


case class Url(url: String)

case class ChecksumObject(checksum: String, `type`: String)

case class DosDataObject(size: Option[Long],
                         checksums: Option[Array[ChecksumObject]],
                         updated: Option[String],
                         urls: Array[Url])

case class DosObject(data_object: DosDataObject)

case class SADataObject(data: Json)

case class MarthaResponse(dos: DosObject, googleServiceAccount: Option[SADataObject])
