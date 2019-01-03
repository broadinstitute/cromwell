package cloud.nio.impl.drs

import java.io.ByteArrayInputStream
import java.nio.channels.{ReadableByteChannel, WritableByteChannel}

import cats.implicits._
import cats.data.{EitherT, Kleisli}
import cats.effect.{Effect, IO}
import cloud.nio.spi.{CloudNioFileList, CloudNioFileProvider, CloudNioRegularFileAttributes}
import com.google.auth.oauth2.{AccessToken, GoogleCredentials, OAuth2Credentials}
import com.google.cloud.ReadChannel
import com.google.cloud.storage.{Storage, StorageOptions}
import com.typesafe.config.Config
import org.apache.http.HttpStatus
import org.apache.http.client.methods.CloseableHttpResponse
import org.apache.http.impl.client.{CloseableHttpClient, HttpClientBuilder}

import scala.util.Try
import scala.concurrent.duration._


class DrsCloudNioFileProvider(config: Config, scheme: String, drsPathResolver: DrsPathResolver, authCredentials: OAuth2Credentials, readInterpreter: (MarthaResponse) => IO[ReadableByteChannel]) extends CloudNioFileProvider {

  private val AccessTokenAcceptableTTL = 1.minute

  private lazy val httpClientBuilder = HttpClientBuilder.create()

  private lazy val marthaUri = config.getString("martha.url")
  private lazy val marthaRequestJsonTemplate = config.getString("martha.request.json-template")

  private val GcsScheme = "gs"


  private def getDrsPath(cloudHost: String, cloudPath: String): String = s"$scheme://$cloudHost/$cloudPath"

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

  private def checkIfPathExistsThroughMartha(drsPath: String): Boolean = {
    val httpClient: CloseableHttpClient = httpClientBuilder.build()

    try {
      val marthaResponse: CloseableHttpResponse = httpClient.execute(drsPathResolver.makeHttpRequestToMartha(drsPath))

      try {
        marthaResponse.getStatusLine.getStatusCode == HttpStatus.SC_OK
      } finally {
        Try(marthaResponse.close())
        ()
      }
    } finally {
      Try(httpClient.close())
      ()
    }
  }


  private def extractUrlForScheme(drsPath: String, urlArray: Array[Url], scheme: String): Either[UrlNotFoundException, String] = {
    val schemeUrlOption = urlArray.find(u => u.url.startsWith(scheme))

    schemeUrlOption match {
      case Some(schemeUrl) => Right(schemeUrl.url)
      case None => Left(UrlNotFoundException(drsPath, scheme))
    }
  }


  private def gcsInputStream(gcsFile: GCSFile): Kleisli[IO, Storage, ReadChannel] =
  Kleisli { storage =>
    IO.delay {
      val blob = storage.get(gcsFile.bucket, gcsFile.file)
      blob.reader()
    }
  }


  private def inputStream(url: String): IO[GCSFile] =  {
    scheme match {
      case GcsScheme => {
        val urlArray = url.replace(s"$GcsScheme://", "").split("/", 2)
        val fileToBeLocalized = urlArray(1)
        val bucket = urlArray(0)
        IO(GCSFile(bucket, fileToBeLocalized))
      }
      case otherScheme => IO.raiseError(throw new UnsupportedOperationException(s"DRS currently doesn't support reading files for $otherScheme."))
    }
  }


  override def existsPath(cloudHost: String, cloudPath: String): Boolean =
    checkIfPathExistsThroughMartha(getDrsPath(cloudHost, cloudPath))


  override def existsPaths(cloudHost: String, cloudPathPrefix: String): Boolean =
    existsPath(cloudHost, cloudPathPrefix)


  override def listObjects(cloudHost: String, cloudPathPrefix: String, markerOption: Option[String]): CloudNioFileList = {
    val exists = existsPath(cloudHost, cloudPathPrefix)
    val list = if (exists) List(cloudPathPrefix) else Nil
    CloudNioFileList(list, None)
  }


  override def copy(sourceCloudHost: String, sourceCloudPath: String, targetCloudHost: String, targetCloudPath: String): Unit =
    throw new UnsupportedOperationException("DRS currently doesn't support copy.")


  override def deleteIfExists(cloudHost: String, cloudPath: String): Boolean =
    throw new UnsupportedOperationException("DRS currently doesn't support delete.")


  override def read(cloudHost: String, cloudPath: String, offset: Long): ReadableByteChannel = {
    val drsPath = getDrsPath(cloudHost,cloudPath)
    val freshAccessToken = getFreshAccessToken(authCredentials)

    val credentials = GoogleCredentials.fromStream(new ByteArrayInputStream(freshAccessToken.getBytes()))
    val storage = StorageOptions.newBuilder().setCredentials(credentials).build().getService





//    val gcsFile: EitherT[String, GCSFile] = inputStream(url)
//    gcsFile.map{
//      file =>
//        gcsInputStream(file).run(storage)
//    } match {
//      case Right(byteChannel) => byteChannel
//      case Left(exception) => throw new RuntimeException(exception)
//
//    }

    val byteChannelIO = for {
      marthaResponse <- drsPathResolver.resolveDrsThroughMartha(drsPath, Option(freshAccessToken))
      byteChannel <- readInterpreter(marthaResponse)

      //Currently, Martha only supports resolving DRS paths to GCS paths
//      url <- IO.fromEither(extractUrlForScheme(drsPath, marthaResponse.dos.data_object.urls, GcsScheme))
//      input <- inputStream(url)
    } yield byteChannel

    byteChannelIO.unsafeRunSync()
  }


  override def write(cloudHost: String, cloudPath: String): WritableByteChannel =
    throw new UnsupportedOperationException("DRS currently doesn't support write.")


  override def fileAttributes(cloudHost: String, cloudPath: String): Option[CloudNioRegularFileAttributes] =
    Option(new DrsCloudNioRegularFileAttributes(config, getDrsPath(cloudHost,cloudPath), drsPathResolver))
}


case class GoogleSANotFoundException(drsPath: String) extends Exception(s"Error finding Google Service Account associated with DRS path $drsPath through Martha.")

case class UrlNotFoundException(drsPath: String, scheme: String) extends Exception(s"DRS was not able to find a $scheme url associated with $drsPath.")

case class GCSFile(bucket: String, file: String)

