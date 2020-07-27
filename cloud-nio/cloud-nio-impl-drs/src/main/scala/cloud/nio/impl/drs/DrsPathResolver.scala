package cloud.nio.impl.drs

import java.time.OffsetDateTime

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
  lazy val marthaResponseV3Decoder: Decoder[MarthaResponseV3] = deriveDecoder

  def marthaV2ToV3ResponseDecoder : Decoder[MarthaResponseV2] = (h: HCursor) => {
    for {
      googleSAOption <- h.downField("googleServiceAccount").as[Option[SADataObject]]
      sizeOption <- h.downField("dos").downField("data_object").downField("size").as[Option[Long]]
      contentTypeOption <- h.downField("dos").downField("data_object").downField("mime_type").as[Option[String]]
      timeCreatedOption <- h.downField("dos").downField("data_object").downField("created").as[Option[String]]
      timeCreatedISO = timeCreatedOption.map(OffsetDateTime.parse(_).toString)
      timeUpdatedOption <- h.downField("dos").downField("data_object").downField("updated").as[Option[String]]
      timeUpdatedISO = timeUpdatedOption.map(OffsetDateTime.parse(_).toString)
      checksumsArray <- h.downField("dos").downField("data_object").downField("checksums").as[List[Map[String, String]]]
      hashesMap = convertToHashesMap(checksumsArray)
      urlsArray <- h.downField("dos").downField("data_object").downField("urls").as[List[Map[String, String]]]
      gcsUrlOption = extractGcsUrl(urlsArray)
      (bucketNameOption, fileNameOption) = getGcsBucketAndName(gcsUrlOption)
    } yield new MarthaResponseV2(
      contentType = contentTypeOption,
      timeCreated = timeCreatedISO,
      timeUpdated = timeUpdatedISO,
      googleServiceAccount = googleSAOption,
      size = sizeOption,
      hashes = Option(hashesMap),
      gsUri = gcsUrlOption,
      bucket = bucketNameOption,
      name = fileNameOption
    )
  }

  private val DrsPathToken = "${drsPath}"


  private def convertToHashesMap(hashesList: List[Map[String, String]]): Map[String, String] = {
    hashesList.flatMap { hashMap =>
      if (hashMap.contains("type") && hashMap.contains("checksum"))
        Map(hashMap("type") -> hashMap("checksum"))
      else
        throw new RuntimeException("Hashes did not contain either `type` or `checksum`.")
    }.toMap
  }


  private def extractGcsUrl(urlList: List[Map[String, String]]): Option[String] = {
    urlList.map { urlMap =>
      val value: Option[String] = urlMap.get("url")
      value match {
        case Some(url) =>
          if (url.startsWith("gs://")) Option(url)
          else None
        case None => None
      }
    }.head
  }

  private def getGcsBucketAndName(gcsUrlOption: Option[String]): (Option[String], Option[String]) = {
    gcsUrlOption match {
      case Some(gcsUrl) =>
        val array = gcsUrl.substring(5).split("/", 1)
        (Option(array(0)), Option(array(1)))
      case None => (None, None)
    }
  }

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

    responseContentIO.flatMap{ responseContent =>
      if (drsConfig.marthaUri.endsWith("martha_v3")) {
        IO.fromEither(decode[MarthaResponseV3](responseContent)(marthaResponseV3Decoder))
      } else {
        IO.fromEither(decode[MarthaResponseV2](responseContent)(marthaV2ToV3ResponseDecoder))
      }
    }.handleErrorWith {
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
                          hashes: Option[Map[String, String]])

class MarthaResponseV2(contentType: Option[String],
                       size: Option[Long],
                       timeCreated: Option[String],
                       timeUpdated: Option[String],
                       bucket: Option[String],
                       name: Option[String],
                       gsUri: Option[String],
                       googleServiceAccount: Option[SADataObject],
                       hashes: Option[Map[String, String]]) extends MarthaResponse(contentType, size, timeCreated, timeUpdated, bucket, name, gsUri, googleServiceAccount, hashes)

class MarthaResponseV3(contentType: Option[String],
                       size: Option[Long],
                       timeCreated: Option[String],
                       timeUpdated: Option[String],
                       bucket: Option[String],
                       name: Option[String],
                       gsUri: Option[String],
                       googleServiceAccount: Option[SADataObject],
                       hashes: Option[Map[String, String]]) extends MarthaResponse(contentType, size, timeCreated, timeUpdated, bucket, name, gsUri, googleServiceAccount, hashes)

