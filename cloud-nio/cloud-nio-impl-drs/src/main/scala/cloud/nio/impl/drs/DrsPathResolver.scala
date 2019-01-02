package cloud.nio.impl.drs

import cats.data.Kleisli
import cats.effect.IO
import com.google.auth.oauth2.{AccessToken, OAuth2Credentials}
import com.typesafe.config.Config
import org.apache.commons.lang3.exception.ExceptionUtils
import org.apache.http.client.methods.{CloseableHttpResponse, HttpPost}
import org.apache.http.entity.{ContentType, StringEntity}
import org.apache.http.impl.client.{CloseableHttpClient, HttpClientBuilder}
import org.apache.http.util.EntityUtils
import org.apache.http.{HttpStatus, StatusLine}
import spray.json._

import scala.util.{Failure, Success, Try}


object DrsPathResolver {
//  private val DrsPathToken = s"$${drsPath}"
//  private lazy val marthaUri = config.getString("martha.url")
//  private lazy val marthaRequestJsonTemplate = config.getString("martha.request.json-template")
//
//  private lazy val httpClientBuilder = HttpClientBuilder.create()

  def makeHttpRequestToMartha(drsPath: String, serviceAccount: Option[String] = None): HttpPost = {
    val postRequest = new HttpPost(marthaUri)
    val requestJson = marthaRequestJsonTemplate.replace(DrsPathToken, drsPath)
    postRequest.setEntity(new StringEntity(requestJson, ContentType.APPLICATION_JSON))
    serviceAccount.foreach(sa => postRequest.setHeader("Authorization", s"Bearer $sa"))
    postRequest
  }


  /***
    * Resolves the DRS path through Martha url provided in the config.
    * Please note, this method makes a synchronous HTTP request to Martha.
    */
  def resolveDrsThroughMartha(drsPath: String, serviceAccount: Option[String] = None): Kleisli[IO, CloseableHttpClient, MarthaResponse] = Kleisli { httpClient =>
    import MarthaResponseJsonSupport._

    def unexpectedExceptionResponse(status: StatusLine) = {
      throw new RuntimeException(s"Unexpected response resolving $drsPath through Martha url $marthaUri. Error: ${status.getStatusCode} ${status.getReasonPhrase}.")
    }


    IO(httpClient.execute(makeHttpRequestToMartha(drsPath, serviceAccount))).bracket { marthaResponse =>
      val marthaResponseEntityOption = Option(marthaResponse.getEntity).map(EntityUtils.toString)
      val responseStatusLine = marthaResponse.getStatusLine

      val responseContentIO = IO.fromEither(marthaResponseEntityOption.filter { _ =>
        responseStatusLine.getStatusCode == HttpStatus.SC_OK
      }.toRight(unexpectedExceptionResponse(responseStatusLine)))

        responseContentIO.flatMap(responseContent =>IO(responseContent.parseJson.convertTo[MarthaResponse])).handleErrorWith{
          e =>
            val a = new RuntimeException(s"Failed to resolve DRS path $drsPath. Error while parsing the response from Martha.")
            a.addSuppressed(e)
            IO.raiseError(a)
        }
    } {httpResponse =>
      IO(httpResponse.close())
    }
  }
}


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
