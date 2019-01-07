package cloud.nio.impl.drs

import cats.effect.{IO, Resource}
import org.apache.http.client.methods.HttpPost
import org.apache.http.entity.{ContentType, StringEntity}
import org.apache.http.impl.client.HttpClientBuilder
import org.apache.http.util.EntityUtils
import org.apache.http.{HttpResponse, HttpStatus, StatusLine}
import spray.json._
//???
import cats.syntax.functor._


case class DrsPathResolver(drsConfig: DrsConfig, httpClientBuilder: HttpClientBuilder) {
  private val DrsPathToken = "$${drsPath}"

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
    import MarthaResponseJsonSupport._

    val marthaResponseEntityOption = Option(httpResponse.getEntity).map(EntityUtils.toString)
    val responseStatusLine = httpResponse.getStatusLine

    val responseContentIO =
      IO.fromEither(
        marthaResponseEntityOption.filter { _ => responseStatusLine.getStatusCode == HttpStatus.SC_OK}
          .toRight(unexpectedExceptionResponse(responseStatusLine))
      )

    responseContentIO.flatMap(responseContent => IO(responseContent.parseJson.convertTo[MarthaResponse])).handleErrorWith {
      e =>
        val errorMsg = new RuntimeException(s"Failed to parse response from Martha into a case class.")
        errorMsg.addSuppressed(e)
        IO.raiseError(errorMsg)
    }
  }


  private def executeMarthaRequest(httpPost: HttpPost): Resource[IO, HttpResponse]= {
//    val a: Resource[IO, CloseableHttpClient] = Resource.fromAutoCloseable(IO(httpClientBuilder.build()))
//    a.flatMap(httpClient =>
//      Resource.fromAutoCloseable(IO(httpClient.execute(httpPost)))
//    )
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


object MarthaResponseJsonSupport extends DefaultJsonProtocol {
  implicit val urlFormat: JsonFormat[Url] = jsonFormat1(Url)
  implicit val checksumFormat: JsonFormat[ChecksumObject] = jsonFormat2(ChecksumObject)
  implicit val dataObject: JsonFormat[DosDataObject] = jsonFormat4(DosDataObject)
  implicit val dosObjectFormat: JsonFormat[DosObject] = jsonFormat1(DosObject)
  implicit val saDataObjectFormat: JsonFormat[SADataObject] = jsonFormat1(SADataObject)
  implicit val marthaResponseFormat: JsonFormat[MarthaResponse] = jsonFormat2(MarthaResponse)
}

case class Url(url: String)

case class ChecksumObject(checksum: String, `type`: String)

case class DosDataObject(size: Option[Long],
                         checksums: Option[Array[ChecksumObject]],
                         updated: Option[String],
                         urls: Array[Url])

case class DosObject(data_object: DosDataObject)

case class SADataObject(data: JsObject)

case class MarthaResponse(dos: DosObject, googleServiceAccount: Option[SADataObject])
