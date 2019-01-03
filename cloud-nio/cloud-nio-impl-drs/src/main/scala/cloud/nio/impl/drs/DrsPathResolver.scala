package cloud.nio.impl.drs

import com.google.auth.oauth2.{AccessToken, OAuth2Credentials}
import com.typesafe.config.Config
import io.circe._
import io.circe.generic.semiauto._
import io.circe.parser.decode
import org.apache.commons.lang3.exception.ExceptionUtils
import org.apache.http.client.methods.{CloseableHttpResponse, HttpPost}
import org.apache.http.entity.{ContentType, StringEntity}
import org.apache.http.impl.client.{CloseableHttpClient, HttpClientBuilder}
import org.apache.http.util.EntityUtils
import org.apache.http.{HttpStatus, StatusLine}

import scala.concurrent.duration._
import scala.util.Try


case class DrsPathResolver(config: Config, authCredentials: OAuth2Credentials) {
  implicit lazy val urlDecoder: Decoder[Url] = deriveDecoder
  implicit lazy val checksumDecoder: Decoder[ChecksumObject] = deriveDecoder
  implicit lazy val dataObjectDecoder: Decoder[DosDataObject] = deriveDecoder
  implicit lazy val dosObjectDecoder: Decoder[DosObject] = deriveDecoder
  implicit lazy val saDataObjectDecoder: Decoder[SADataObject] = deriveDecoder
  implicit lazy val marthaResponseDecoder: Decoder[MarthaResponse] = deriveDecoder

  private val AccessTokenAcceptableTTL = 1.minute

  private val DrsPathToken = s"$${drsPath}"
  private lazy val marthaUri = config.getString("martha.url")
  private lazy val marthaRequestJsonTemplate = config.getString("martha.request.json-template")

  private lazy val httpClientBuilder = HttpClientBuilder.create()


  private def getFreshAccessToken(credential: OAuth2Credentials): String = {
    def accessTokenTTLIsAcceptable(accessToken: AccessToken): Boolean = {
      (accessToken.getExpirationTime.getTime - System.currentTimeMillis()).millis.gteq(AccessTokenAcceptableTTL)
    }

    Option(credential.getAccessToken) match {
      case Some(accessToken) if accessTokenTTLIsAcceptable(accessToken) => accessToken.getTokenValue
      case _ =>
        credential.refresh()
        credential.getAccessToken.getTokenValue
    }
  }


  def makeHttpRequestToMartha(drsPath: String, httpClient: CloseableHttpClient, needServiceAccount: Boolean): CloseableHttpResponse = {
    val postRequest = new HttpPost(marthaUri)
    val requestJson = marthaRequestJsonTemplate.replace(DrsPathToken, drsPath)
    postRequest.setEntity(new StringEntity(requestJson, ContentType.APPLICATION_JSON))

    if (needServiceAccount) {
      val accessToken = getFreshAccessToken(authCredentials)
      postRequest.setHeader("Authorization", s"Bearer $accessToken")
    }

    httpClient.execute(postRequest)
  }


  /***
    * Resolves the DRS path through Martha url provided in the config.
    * Please note, this method makes a synchronous HTTP request to Martha.
    */
  def resolveDrsThroughMartha(drsPath: String, needServiceAccount: Boolean = false): MarthaResponse = {

    def unexpectedExceptionResponse(status: StatusLine) = {
      throw new RuntimeException(s"Unexpected response resolving $drsPath through Martha url $marthaUri. Error: ${status.getStatusCode} ${status.getReasonPhrase}.")
    }

    val httpClient: CloseableHttpClient = httpClientBuilder.build()

    try {
      val marthaResponse: CloseableHttpResponse = makeHttpRequestToMartha(drsPath, httpClient, needServiceAccount)
      val marthaResponseEntityOption = Option(marthaResponse.getEntity).map(EntityUtils.toString)

      try {
        val responseStatusLine = marthaResponse.getStatusLine
        val responseContent = if (responseStatusLine.getStatusCode == HttpStatus.SC_OK) {
          marthaResponseEntityOption.getOrElse(unexpectedExceptionResponse(responseStatusLine))
        } else {
          unexpectedExceptionResponse(responseStatusLine)
        }

        decode[MarthaResponse](responseContent) match {
          case Right(marthaObj) => marthaObj
          case Left(e) => throw new RuntimeException(s"Failed to resolve DRS path $drsPath. Error while parsing the response from Martha. Error: ${ExceptionUtils.getMessage(e)}")
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


case class Url(url: String)

case class ChecksumObject(checksum: String, `type`: String)

case class DosDataObject(size: Option[Long],
                         checksums: Option[Array[ChecksumObject]],
                         updated: Option[String],
                         urls: Array[Url])

case class DosObject(data_object: DosDataObject)

case class SADataObject(data: JsonObject)

case class MarthaResponse(dos: DosObject, googleServiceAccount: Option[SADataObject])
