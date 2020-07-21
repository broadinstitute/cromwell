package cloud.nio.impl.drs

import cats.effect.{IO, Resource}
import cats.instances.option._
import cats.instances.string._
import com.google.auth.oauth2.{AccessToken, OAuth2Credentials}
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

import scala.concurrent.duration._

case class DrsPathResolver(drsConfig: DrsConfig,
                           httpClientBuilder: HttpClientBuilder,
                           accessTokenAcceptableTTL: Duration,
                           authCredentials: OAuth2Credentials) {

  implicit lazy val saDataObjectDecoder: Decoder[SADataObject] = deriveDecoder
  implicit lazy val marthaResponseDecoder: Decoder[MarthaResponse] = deriveDecoder

  private val DrsPathToken = "${drsPath}"

  //Based on method from GcrRegistry
  private def getAccessToken: String = {
    def accessTokenTTLIsAcceptable(accessToken: AccessToken): Boolean = {
      (accessToken.getExpirationTime.getTime - System.currentTimeMillis()).millis.gteq(accessTokenAcceptableTTL)
    }

    Option(authCredentials.getAccessToken) match {
      case Some(accessToken) if accessTokenTTLIsAcceptable(accessToken) => accessToken.getTokenValue
      case _ =>
        authCredentials.refresh()
        authCredentials.getAccessToken.getTokenValue
    }
  }


  private def makeHttpRequestToMartha(drsPath: String): HttpPost = {
    val postRequest = new HttpPost(drsConfig.marthaUri)
    val requestJson = drsConfig.marthaRequestJsonTemplate.replace(DrsPathToken, drsPath)
    postRequest.setEntity(new StringEntity(requestJson, ContentType.APPLICATION_JSON))
    postRequest.setHeader("Authorization", s"Bearer $getAccessToken")
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


  def rawMarthaResponse(drsPath: String): Resource[IO, HttpResponse] = {
    val httpPost = makeHttpRequestToMartha(drsPath)
    executeMarthaRequest(httpPost)
  }

  /** *
    * Resolves the DRS path through Martha url provided in the config.
    * Please note, this method returns an IO that would make a synchronous HTTP request to Martha when run.
    */
  def resolveDrsThroughMartha(drsPath: String): IO[MarthaResponse] = {
    rawMarthaResponse(drsPath).use(httpResponseToMarthaResponse)
  }
}

case class DrsConfig(marthaUri: String, marthaRequestJsonTemplate: String)

case class SADataObject(data: Json)

case class MarthaResponse(contentType: Option[String],
                          size: Option[Long],
                          timeCreated: Option[String],
                          timeUpdated: Option[String],
                          bucket: Option[String],
                          name: Option[String],
                          gsUri: Option[String],
                          googleServiceAccount: Option[SADataObject],
                          hashes: Option[JsonObject])
