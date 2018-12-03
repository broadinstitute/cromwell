package cloud.nio.impl.drs

import com.typesafe.config.Config
import org.apache.http.{HttpStatus, StatusLine}
import org.apache.http.client.methods.{CloseableHttpResponse, HttpPost}
import org.apache.http.entity.{ContentType, StringEntity}
import org.apache.http.impl.client.{CloseableHttpClient, HttpClientBuilder}
import org.apache.commons.lang3.exception.ExceptionUtils
import org.apache.http.util.EntityUtils
import spray.json._

import scala.util.{Failure, Success, Try}


case class DrsPathResolver(config: Config) {
  private val DrsPathToken = s"$${drsPath}"
  private lazy val marthaUri = config.getString("martha.url")
  private lazy val marthaRequestJsonTemplate = config.getString("martha.request.json-template")

  //RENAME THIS??
  def contactMartha(drsPath: String, httpClient: CloseableHttpClient): CloseableHttpResponse = {
    val postRequest = new HttpPost(marthaUri)
    val requestJson = marthaRequestJsonTemplate.replace(DrsPathToken, drsPath)
    postRequest.setEntity(new StringEntity(requestJson, ContentType.APPLICATION_JSON))
    httpClient.execute(postRequest)
  }


  def resolveDrsThroughMartha(drsPath: String): MarthaResponse = {
    import MarthaResponseJsonSupport._

    def unexpectedExceptionResponse(status: StatusLine) = {
      throw new RuntimeException(s"Unexpected response resolving $drsPath through Martha url $marthaUri. Error: ${status.getStatusCode} ${status.getReasonPhrase}.")
    }

    val httpClient: CloseableHttpClient = HttpClientBuilder.create().build()

    try {
      val marthaResponse: CloseableHttpResponse = contactMartha(drsPath, httpClient)
      val marthaResponseEntityOption = Option(marthaResponse.getEntity).map(EntityUtils.toString)

      try {
        val responseStatusLine = marthaResponse.getStatusLine
        val responseContent = if (responseStatusLine.getStatusCode == HttpStatus.SC_OK) {
          marthaResponseEntityOption.getOrElse(unexpectedExceptionResponse(responseStatusLine))
        } else {
          unexpectedExceptionResponse(responseStatusLine)
        }

        Try(responseContent.parseJson.convertTo[MarthaResponse]) match {
          case Success(marthaObj) =>
            marthaObj
          case Failure(e) =>
            throw new RuntimeException(s"Failed to resolve DRS path $drsPath. Error while parsing the response from Martha. Error: ${ExceptionUtils.getMessage(e)}")
        }
      } finally {
        Try(marthaResponse.close())
        ()
      }
    } finally {
      Try(httpClient.close())
      ()
    }
  }
}


object MarthaResponseJsonSupport extends DefaultJsonProtocol {
  implicit val urlFormat: JsonFormat[Url] = jsonFormat1(Url)
  implicit val checksumFormat: JsonFormat[ChecksumObject] = jsonFormat2(ChecksumObject)
  implicit val dataObject: JsonFormat[DosDataObject] = jsonFormat4(DosDataObject)
  implicit val dosObjectFormat: JsonFormat[DosObject] = jsonFormat1(DosObject)
  implicit val marthaResponseFormat: JsonFormat[MarthaResponse] = jsonFormat1(MarthaResponse)
}

case class Url(url: String)

case class ChecksumObject(checksum: String, `type`: String)

case class DosDataObject(size: Option[Long],
                         checksums: Option[Array[ChecksumObject]],
                         updated: Option[String],
                         urls: Array[Url])

case class DosObject(data_object: DosDataObject)

case class MarthaResponse(dos: DosObject)
