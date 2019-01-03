#!/usr/bin/env amm

interp.repositories() ++= Seq(
  coursier.maven.MavenRepository("https://oss.sonatype.org/content/repositories/snapshots/"))

@
import $ivy.`org.http4s::http4s-dsl:0.19.0-SNAPSHOT`
import $ivy.`org.http4s::http4s-blaze-server:0.19.0-SNAPSHOT`
import $ivy.`org.http4s::http4s-blaze-client:0.19.0-SNAPSHOT`
import $ivy.`com.google.cloud:google-cloud-storage:1.35.0`
import $ivy.`com.google.oauth-client:google-oauth-client:1.23.0`
import $ivy.`com.google.auth:google-auth-library-credentials:0.9.1`
import $ivy.`org.http4s::http4s-circe:0.19.0-SNAPSHOT`
import $ivy.`io.circe::circe-literal:0.7.0`
import $ivy.`io.spray::spray-json:1.3.4`

import com.google.auth.oauth2.GoogleCredentials
import com.google.cloud.storage.Storage
import com.google.cloud.storage.StorageOptions
import com.google.common.collect.Lists
import com.google.cloud.ReadChannel
import java.io.ByteArrayInputStream
import java.io.FileInputStream
import java.io.FileOutputStream
import java.nio.file.Files
import java.nio.file.Paths
import org.http4s.client.blaze._
import cats.effect._
import org.http4s.dsl.io._
import org.http4s._
import org.http4s.circe._
import io.circe.literal._
import io.circe.syntax._
import spray.json._
import spray.json.DefaultJsonProtocol._
import scala.util.{Try, Success, Failure}

object MarthaResponseJsonSupport extends DefaultJsonProtocol {
  implicit val urlFormat: JsonFormat[Url] = jsonFormat1(Url)
  implicit val dataObject: JsonFormat[DosDataObject] = jsonFormat1(DosDataObject)
  implicit val dosObjectFormat: JsonFormat[DosObject] = jsonFormat1(DosObject)
  implicit val googleServiceAccountFormat: JsonFormat[GoogleServiceAccount] = jsonFormat1(GoogleServiceAccount)
  implicit val marthaResponseFormat: JsonFormat[MarthaResponse] = jsonFormat2(MarthaResponse)
}

case class Url(url: String)

case class DosDataObject(urls: Array[Url])

case class DosObject(data_object: DosDataObject)

case class GoogleServiceAccount(data: JsObject)

case class MarthaResponse(dos: DosObject, googleServiceAccount: GoogleServiceAccount)


@main
def dosUrlResolver(dosUrl: String, downloadLoc: String) : Unit = {
  val dosResloverObj = for {
    marthaUrl <- Uri.fromString(sys.env("MARTHA_URL")).toTry
    marthaResObj <- resolveDosThroughMartha(dosUrl, marthaUrl)
    gcsUrl <- extractFirstGcsUrl(marthaResObj.dos.data_object.urls)
    _ <- downloadFileFromGcs(gcsUrl, marthaResObj.googleServiceAccount.data.toString, downloadLoc)
  } yield()

  dosResloverObj match {
    case Success(_) =>
    case Failure(e) => {
      Console.err.println(s"Error: $e")
      e.printStackTrace(Console.err)
      System.exit(1)
    }
  }
}


def resolveDosThroughMartha(dosUrl: String, marthaUrl: Uri) : Try[MarthaResponse] = {
  import MarthaResponseJsonSupport._

  val requestBody = json"""{"url":$dosUrl}"""

  val credentials = GoogleCredentials.getApplicationDefault()
  val accessToken = credentials.refreshAccessToken().getTokenValue()

  val marthaResponseIo: IO[MarthaResponse] = for {
    httpClient <- Http1Client[IO]()
    postRequest <- Request[IO](method = Method.POST,
                               uri = marthaUrl,
                               headers = Headers(Header("Authorization", s"bearer $accessToken")))
                              .withBody(requestBody)
    httpResponse <- httpClient.expect[String](postRequest)
    marthaResObj = httpResponse.parseJson.convertTo[MarthaResponse]
  } yield marthaResObj

  Try(marthaResponseIo.unsafeRunSync())
}


def extractFirstGcsUrl(urlArray: Array[Url]): Try[String] = {
  val urlObjOption = urlArray.find(urlObj => urlObj.url.startsWith("gs://"))

  urlObjOption match {
    case Some(urlObj) => Success(urlObj.url)
    case None => Failure(new Exception("No resolved url starting with 'gs://' found from Martha response!"))
  }
}


def downloadFileFromGcs(gcsUrl: String, serviceAccount: String, downloadLoc: String) : Try[Unit] = {
  val gcsUrlArray = gcsUrl.replace("gs://", "").split("/", 2)
  val fileToBeLocalized = gcsUrlArray(1)
  val gcsBucket = gcsUrlArray(0)

  Try {
    val credentials = GoogleCredentials.fromStream(new ByteArrayInputStream(serviceAccount.getBytes()))
      .createScoped(Lists.newArrayList("https://www.googleapis.com/auth/cloud-platform"))
    val storage = StorageOptions.newBuilder().setCredentials(credentials).build().getService()
    val blob = storage.get(gcsBucket, fileToBeLocalized)
    val readChannel = blob.reader()
    Files.createDirectories(Paths.get(downloadLoc).getParent)
    val fileOuputStream = new FileOutputStream(downloadLoc)
    fileOuputStream.getChannel().transferFrom(readChannel, 0, Long.MaxValue)
    fileOuputStream.close()
  }
}
