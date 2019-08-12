package drs.localizer

import java.io.ByteArrayInputStream
import java.nio.file.{Files, Paths}

import cats.effect.{ExitCode, IO, IOApp}
import io.circe.parser.decode
import com.google.auth.oauth2.GoogleCredentials
import com.google.cloud.storage.Storage.BlobGetOption
import com.google.cloud.storage.{Blob, StorageException, StorageOptions}
import com.softwaremill.sttp._
import com.softwaremill.sttp.circe._
import drs.localizer.MarthaResponseJsonSupport._
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.concurrent.duration._

object DrsLocalizerMain extends IOApp {

  implicit val httpBackendConnection = HttpURLConnectionBackend(options = SttpBackendOptions.connectionTimeout(5.minutes))

  val logger = LoggerFactory.getLogger("DrsLocalizerLogger")

  val GcsScheme = "gs://"
  val RequesterPaysErrorMsg = "Bucket is requester pays bucket but no user project provided."

  val CloudPlatformAuthScope = "https://www.googleapis.com/auth/cloud-platform"
  val UserInfoEmailScope = "https://www.googleapis.com/auth/userinfo.email"
  val UserInfoProfileScope = "https://www.googleapis.com/auth/userinfo.profile"
  val UserInfoScopes = List(UserInfoEmailScope, UserInfoProfileScope)


  def resolveDrsThroughMartha(drsUrl: String, marthaUrl: Uri): IO[MarthaResponse] = {
    val requestBody = raw"""{"url":"$drsUrl"}"""
    val marthaErrorMsg = s"Something went wrong while trying to resolve $drsUrl through Martha $marthaUrl."

    val scopedCredentials = GoogleCredentials.getApplicationDefault().createScoped(UserInfoScopes.asJava)
    val accessToken = scopedCredentials.refreshAccessToken().getTokenValue

    val request = sttp
      .headers(Map(HeaderNames.ContentType -> MediaTypes.Json, "Authorization" -> s"Bearer $accessToken"))
      .body(requestBody)
      .post(marthaUrl)
      .response(asJson[MarthaResponse])
      .readTimeout(5.minutes)

    val response = request.send()

    response.body match {
      // non-2xx status code (https://sttp.readthedocs.io/en/latest/responses/basics.html)
      case Left(error) => {
        // Extract response from Sam. Martha usually responds with 502, but the response body does contains
        // response from Sam. Extract that for a more helpful error message
        val samErrorResponse = decode[MarthaErrorResponse](error) match {
          case Left(samError) => s"Unable to parse response body from Sam. Error: $samError"
          case Right(samResponse) => s"Response from Sam: Status- ${samResponse.status}. Error message- \n${samResponse.response.text}"
        }
        IO.raiseError(new Exception(s"$marthaErrorMsg Expected 200 but got ${response.code} from Martha. $samErrorResponse"))
      }
      case Right(marthaResponseEither) => {
        logger.info("Received successful response from Martha")
        marthaResponseEither match {
          case Left(deserializationError) => IO.raiseError(new Exception(s"$marthaErrorMsg Failed to parse Martha response. Deserialization error: ${deserializationError.message}"))
          case Right(marthaResponse) => IO(marthaResponse)
        }
      }
    }
  }


  def extractFirstGcsUrl(urlArray: Array[Url]): IO[String] = {
    val urlOption = urlArray.find(urlObj => urlObj.url.startsWith(GcsScheme))

    urlOption match {
      case Some(url) => IO(url.url)
      case None => IO.raiseError(new Exception("No resolved url starting with 'gs://' found from Martha response!"))
    }
  }


  def downloadFileFromGcs(gcsUrl: String,
                          serviceAccountJsonOption: Option[String],
                          downloadLoc: String,
                          requesterPaysProjectIdOption: Option[String]) : IO[Unit] = {
    logger.info(s"Requester Pays project ID is $requesterPaysProjectIdOption")

    val gcsUrlArray = gcsUrl.replace(GcsScheme, "").split("/", 2)
    val fileToBeLocalized = gcsUrlArray(1)
    val gcsBucket = gcsUrlArray(0)

    val unscopedCredentials = serviceAccountJsonOption match {
      case None => GoogleCredentials.getApplicationDefault()
      case Some(serviceAccountJson) => GoogleCredentials.fromStream(new ByteArrayInputStream(serviceAccountJson.getBytes()))
    }

    val credentials = unscopedCredentials.createScoped(List(CloudPlatformAuthScope).asJava)

    val storage = StorageOptions.newBuilder().setCredentials(credentials).build().getService
    Files.createDirectories(Paths.get(downloadLoc).getParent)

    logger.info(s"Attempting to download $gcsUrl")

    IO.delay {
      val blob = storage.get(gcsBucket, fileToBeLocalized)
      blob.downloadTo(Paths.get(downloadLoc))
      logger.info(s"Download complete without using Requester Pays")
    } handleErrorWith { throwable =>
      (requesterPaysProjectIdOption, throwable) match {
        case (Some(requesterPaysProjectId), storageException: StorageException) if storageException.getMessage == RequesterPaysErrorMsg => {
          logger.info(s"Received 'Bucket is requester pays' error. Attempting again using Requester Pays billing project")
          IO.delay {
            val blob = storage.get(gcsBucket, fileToBeLocalized, BlobGetOption.userProject(requesterPaysProjectId))
            blob.downloadTo(Paths.get(downloadLoc), Blob.BlobSourceOption.userProject(requesterPaysProjectId))
          }
        }
        case _ => IO.raiseError(throwable)
      }
    }
  }


  //args(0) - dos url
  //args(1) - download location
  //args(2) - martha url for now
  override def run(args: List[String]): IO[ExitCode] = {
    val existStateIO = for {
      marthaUri <- IO(uri"${args(2)}")
      marthaResponse <- resolveDrsThroughMartha(args.head, marthaUri)
      _ = httpBackendConnection.close()
      gcsUrl <- extractFirstGcsUrl(marthaResponse.dos.data_object.urls)
      exitState <- downloadFileFromGcs(gcsUrl, marthaResponse.googleServiceAccount.map(_.data.toString), args(1), None) //TODO: Saloni- this should not be a NONE!!
    } yield exitState

    existStateIO.map(_ => ExitCode.Success)
  }
}
