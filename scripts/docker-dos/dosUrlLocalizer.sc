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
import java.lang.Long
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


case class MarthaResponse(dos: String, sa: String)

object MarthaResponseJsonSupport extends DefaultJsonProtocol {
  implicit val responseFormat = jsonFormat2(MarthaResponse)
}


@main
def dosUrlResolver(dosUrl: String, downloadLoc: String) : Unit = {
  val dosResloverObj = for {
    marthaResObj <- resolveDosThroughMartha(dosUrl)
    _ <- downloadFileFromGcs(marthaResObj.dos, marthaResObj.sa, downloadLoc)
  } yield()

  dosResloverObj match {
    case Success(_) =>
    case Failure(e) => {
      println("Error:")
      println(e.printStackTrace)
      System.exit(1)
    }
  }
}


def resolveDosThroughMartha(dosUrl: String) : Try[MarthaResponse] = {
  import MarthaResponseJsonSupport._

  // if using fake Martha request, insert actual service account json here
  // if using actual Martha remove this variable entirely
  val serviceAccount = raw"""{}"""

  val requestBody = json"""{"dosUrl":$dosUrl}"""

  Try {
    val httpClient = Http1Client[IO]().unsafeRunSync
    //request to fake Martha
    val postRequest = Request[IO](
      method = Method.POST,
      uri = Uri.uri("https://us-central1-broad-dsde-cromwell-dev.cloudfunctions.net/helloWorld"))
      .withBody(requestBody)
    val httpResponse = httpClient.expect[String](postRequest).unsafeRunSync
    val marthaResObj = httpResponse.parseJson.convertTo[MarthaResponse]

    //reconstructing the marthaResObj to replace sa details
    MarthaResponse(dos = marthaResObj.dos, sa = serviceAccount)
  }
}


def downloadFileFromGcs(gcsUrl: String, serviceAccount: String, downloadLoc: String) : Try[Unit] = {
  val gcsUrlArray = gcsUrl.split("/")
  val fileToBeLocalized = gcsUrlArray(3)
  val gcsBucket = gcsUrlArray(2)

  for {
    credentials <- Try(GoogleCredentials.fromStream(new ByteArrayInputStream(serviceAccount.getBytes()))
      .createScoped(Lists.newArrayList("https://www.googleapis.com/auth/cloud-platform")))
    storage = StorageOptions.newBuilder().setCredentials(credentials).build().getService()
    blob <- Try(storage.get(gcsBucket, fileToBeLocalized))
    readChannel <- Try(blob.reader())
    _ <- Try(Files.createDirectories(Paths.get(downloadLoc).getParent))
    fileOuputStream = new FileOutputStream(downloadLoc)
    _ <- Try(fileOuputStream.getChannel().transferFrom(readChannel, 0, Long.MAX_VALUE))
    _ <- Try(fileOuputStream.close())
  } yield()
}
