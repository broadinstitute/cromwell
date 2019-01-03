package cloud.nio.impl.drs

import cats.effect.{IO, Resource}
import org.apache.http.client.methods.HttpPost
import org.apache.http.entity.{ContentType, StringEntity}
import org.apache.http.impl.client.CloseableHttpClient
import org.apache.http.util.EntityUtils
import org.apache.http.{HttpStatus, StatusLine}
import spray.json._


case class DrsPathResolver(drsConfig: DrsConfig, closeableHttpClient: CloseableHttpClient) {
  private val DrsPathToken = "$${drsPath}"

//  private lazy val httpClientBuilder = HttpClientBuilder.create()

  def makeHttpRequestToMartha(drsPath: String, serviceAccount: Option[String] = None): IO[HttpPost] = {
    val postRequest = new HttpPost(drsConfig.marthaUri)
    val requestJson = drsConfig.marthaRequestJsonTemplate.replace(DrsPathToken, drsPath)
    postRequest.setEntity(new StringEntity(requestJson, ContentType.APPLICATION_JSON))
    serviceAccount.foreach(sa => postRequest.setHeader("Authorization", s"Bearer $sa"))
    IO(postRequest)
  }


  def executeMarthaRequest(httpPost: HttpPost): IO[MarthaResponse] = {
    import MarthaResponseJsonSupport._

    def unexpectedExceptionResponse(status: StatusLine) = {
      new RuntimeException(s"Unexpected response resolving DRS path through Martha url ${drsConfig.marthaUri}. Error: ${status.getStatusCode} ${status.getReasonPhrase}.")
    }

    Resource.fromAutoCloseable(IO(closeableHttpClient.execute(httpPost))).use { marthaResponse =>
      val marthaResponseEntityOption = Option(marthaResponse.getEntity).map(EntityUtils.toString)
      val responseStatusLine = marthaResponse.getStatusLine

      val responseContentIO = IO.fromEither(marthaResponseEntityOption.filter { _ =>
        responseStatusLine.getStatusCode == HttpStatus.SC_OK
      }.toRight(unexpectedExceptionResponse(responseStatusLine)))

      responseContentIO.flatMap(responseContent =>IO(responseContent.parseJson.convertTo[MarthaResponse])).handleErrorWith{
        e =>
          val a = new RuntimeException(s"Falied to parse JSON into Martha response.")
          a.addSuppressed(e)
          IO.raiseError(a)
      }
    }
  }

  /***
    * Resolves the DRS path through Martha url provided in the config.
    * Please note, this method makes a synchronous HTTP request to Martha.
    */
  def resolveDrsThroughMartha(drsPath: String, serviceAccount: Option[String] = None): IO[MarthaResponse] = {
    for {
      httpPost <- makeHttpRequestToMartha(drsPath, serviceAccount)
      marthaResponse <- executeMarthaRequest(httpPost)
    } yield marthaResponse
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
