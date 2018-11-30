package cloud.nio.impl.drs

import com.typesafe.config.Config
import org.apache.http.HttpStatus
import org.apache.http.client.methods.{CloseableHttpResponse, HttpPost}
import org.apache.http.entity.{ContentType, StringEntity}
import org.apache.http.impl.client.{CloseableHttpClient, HttpClientBuilder}
import org.apache.http.util.EntityUtils
import spray.json._

import scala.util.Try


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

    def unexpectedExceptionResponse = {
      throw new RuntimeException(s"Unexpected response looking up $drsPath from $marthaUri.")
    }

    val httpClient: CloseableHttpClient = HttpClientBuilder.create().build()

    try {
      val marthaResponse: CloseableHttpResponse = contactMartha(drsPath, httpClient)
      val marthaResponseEntityOption = Option(marthaResponse.getEntity).map(EntityUtils.toString)

      try {
        val responseContent = if (marthaResponse.getStatusLine.getStatusCode == HttpStatus.SC_OK) {
          marthaResponseEntityOption.getOrElse(unexpectedExceptionResponse)
        } else {
          unexpectedExceptionResponse
        }

        responseContent.parseJson.convertTo[MarthaResponse]
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
  //  implicit val googleServiceAccountFormat: JsonFormat[GoogleServiceAccount] = jsonFormat1(GoogleServiceAccount)
  implicit val marthaResponseFormat: JsonFormat[MarthaResponse] = jsonFormat2(MarthaResponse)
}

case class Url(url: String)

case class ChecksumObject(checksum: String, `type`: String)

case class DosDataObject(size: Long,
                         checksums: Array[ChecksumObject],
                         updated: String,
                         urls: Array[Url])

case class DosObject(data_object: DosDataObject)

//DO NOT MERGE IF BELOW CODE IS COMMENTED OUT
//The Martha response from prod/dev url does contain `data` key inside googleServiceAccount,
//but currently for testing purposes since mock martha doesn't have `data` key it is commented out
//case class GoogleServiceAccount(data: JsObject)
//case class MarthaResponse(dos: DosObject, googleServiceAccount: GoogleServiceAccount)


case class MarthaResponse(dos: DosObject, googleServiceAccount: JsObject)
