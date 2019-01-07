package cloud.nio.impl.drs

import java.nio.channels.{ReadableByteChannel, WritableByteChannel}

import cats.effect.IO
import cloud.nio.spi.{CloudNioFileList, CloudNioFileProvider, CloudNioRegularFileAttributes}
import com.google.auth.oauth2.{AccessToken, OAuth2Credentials}
import com.typesafe.config.Config
import org.apache.http.HttpStatus
import org.apache.http.impl.client.HttpClientBuilder

import scala.concurrent.duration._


class DrsCloudNioFileProvider(config: Config,
                              scheme: String,
                              drsPathResolver: DrsPathResolver,
                              authCredentials: OAuth2Credentials,
                              httpClientBuilder: HttpClientBuilder,
                              drsReadInterpreter: (String, MarthaResponse) => IO[ReadableByteChannel]) extends CloudNioFileProvider {

  private val AccessTokenAcceptableTTL = 1.minute

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
    val a = drsPathResolver.rawMarthaResponse(drsPath).use { marthaResponse =>
      IO.fromEither(Option(marthaResponse.getStatusLine).toRight(new RuntimeException(s"Status line was null for martha response $marthaResponse.")))
    }.map(_.getStatusCode == HttpStatus.SC_OK)

    a.unsafeRunSync()
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

    val byteChannelIO = for {
      marthaResponse <- drsPathResolver.resolveDrsThroughMartha(drsPath, Option(freshAccessToken))
      byteChannel <- drsReadInterpreter(getDrsPath(cloudHost, cloudPath), marthaResponse)
    } yield byteChannel

    byteChannelIO.unsafeRunSync()
  }


  override def write(cloudHost: String, cloudPath: String): WritableByteChannel =
    throw new UnsupportedOperationException("DRS currently doesn't support write.")


  override def fileAttributes(cloudHost: String, cloudPath: String): Option[CloudNioRegularFileAttributes] =
    Option(new DrsCloudNioRegularFileAttributes(config, getDrsPath(cloudHost,cloudPath), drsPathResolver))
}



case class GcsFilePath(bucket: String, file: String)

