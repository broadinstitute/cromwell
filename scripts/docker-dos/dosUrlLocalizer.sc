#!/usr/bin/env amm

interp.repositories() ++= Seq(
  coursier.maven.MavenRepository("https://oss.sonatype.org/content/repositories/snapshots/"))

@
import $ivy.`org.http4s::http4s-dsl:0.18.17`
import $ivy.`org.http4s::http4s-blaze-server:0.18.17`
import $ivy.`org.http4s::http4s-blaze-client:0.18.17`
import $ivy.`com.google.cloud:google-cloud-storage:1.35.0`
import $ivy.`com.google.oauth-client:google-oauth-client:1.23.0`
import $ivy.`com.google.auth:google-auth-library-credentials:0.9.1`
import $ivy.`org.http4s::http4s-circe:0.18.17`
import $ivy.`io.circe::circe-literal:0.7.0`
import $ivy.`io.spray::spray-json:1.3.4`

import com.google.auth.oauth2.GoogleCredentials
import com.google.cloud.storage.Storage
import com.google.cloud.storage.StorageException
import com.google.cloud.storage.StorageOptions
import com.google.cloud.storage.Storage.BlobGetOption
import com.google.common.collect.Lists
import com.google.cloud.ReadChannel
import com.google.cloud.storage.Blob
import java.io.ByteArrayInputStream
import java.io.FileInputStream
import java.io.FileOutputStream
import java.nio.file.Files
import java.nio.file.Paths
import java.util.ArrayList
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
import scala.concurrent.duration._
import org.http4s.client._


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

case class MarthaResponse(dos: DosObject, googleServiceAccount: Option[GoogleServiceAccount])


@main
def dosUrlResolver(dosUrl: String, downloadLoc: String) : Unit = {
  val dosResloverObj = for {
    marthaUrl <- Uri.fromString(sys.env("MARTHA_URL")).toTry
    requesterPaysProjectIdOption <- Try(sys.env.get("REQUESTER_PAYS_PROJECT_ID"))
    marthaResObj <- resolveDosThroughMartha(dosUrl, marthaUrl)
    gcsUrl <- extractFirstGcsUrl(marthaResObj.dos.data_object.urls)
    _ <- downloadFileFromGcs(
      gcsUrl,
      marthaResObj.googleServiceAccount.map(_.data.toString),
      downloadLoc,
      requesterPaysProjectIdOption
    )
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


/*
* Note: The printlns introduced are temporary and would help with debugging problems while running this script.
* They should be replaced by logging messages in ticket: https://broadworkbench.atlassian.net/browse/BA-5821.
*/
def resolveDosThroughMartha(dosUrl: String, marthaUrl: Uri) : Try[MarthaResponse] = {
  import MarthaResponseJsonSupport._

  println(s"Resolving $dosUrl to martha $marthaUrl")

  val requestBody = json"""{"url":$dosUrl}"""
  val scopes = Lists.newArrayList("https://www.googleapis.com/auth/userinfo.email",
    "https://www.googleapis.com/auth/userinfo.profile")

  val credentials = GoogleCredentials.getApplicationDefault()
  val scopedCredentials = credentials.createScoped(scopes)
  val accessToken = scopedCredentials.refreshAccessToken().getTokenValue()

  val longTimeoutConfig =
    BlazeClientConfig
      .defaultConfig
      .copy(idleTimeout = 5.minutes,
        responseHeaderTimeout = 5.minutes,
        requestTimeout = 5.minutes)

  val marthaResponseIo: IO[MarthaResponse] = for {
    httpClient <- Http1Client[IO](config = longTimeoutConfig)
    postRequest <- Request[IO](method = Method.POST,
      uri = marthaUrl,
      headers = Headers(Header("Authorization", s"bearer $accessToken")))
      .withBody(requestBody)
    httpResponse <- httpClient.expect[String](postRequest)
    _ = println("Received successful response from Martha")
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


def downloadFileFromGcs(gcsUrl: String,
                        serviceAccountJsonOption: Option[String],
                        downloadLoc: String,
                        requesterPaysProjectIdOption: Option[String]) : Try[Unit] = {
  println(s"Requester Pays project ID is $requesterPaysProjectIdOption")

  val gcsUrlArray = gcsUrl.replace("gs://", "").split("/", 2)
  val fileToBeLocalized = gcsUrlArray(1)
  val gcsBucket = gcsUrlArray(0)

  Try {
    val unscopedCredentials = serviceAccountJsonOption match {
      case None => GoogleCredentials.getApplicationDefault()
      case Some(serviceAccountJson) =>
        GoogleCredentials.fromStream(new ByteArrayInputStream(serviceAccountJson.getBytes()))
    }
    val credentials =
      unscopedCredentials.createScoped(Lists.newArrayList("https://www.googleapis.com/auth/cloud-platform"))
    val storage = StorageOptions.newBuilder().setCredentials(credentials).build().getService()
    Files.createDirectories(Paths.get(downloadLoc).getParent)

    println(s"Attempting to download $gcsUrl")

    try {
      val blob = storage.get(gcsBucket, fileToBeLocalized)
      blob.downloadTo(Paths.get(downloadLoc))
      println(s"Download complete without using Requester Pays $downloadLoc ")
    } catch {
      case storageException: StorageException
        if storageException.getMessage == "Bucket is requester pays bucket but no user project provided." &&
          requesterPaysProjectIdOption.nonEmpty =>
        println(s"Received StorageException with message stating 'Bucket is requester pays bucket but no user project provided'. " +
          s"Attempting to download with billing project ID")

        requesterPaysProjectIdOption foreach { requesterPaysProjectId =>
          val blob = storage.get(gcsBucket, fileToBeLocalized, BlobGetOption.userProject(requesterPaysProjectId))
          blob.downloadTo(Paths.get(downloadLoc), Blob.BlobSourceOption.userProject(requesterPaysProjectId))
        }

        println(s"Download complete with Requester Pays $downloadLoc")
    }
  }
}
