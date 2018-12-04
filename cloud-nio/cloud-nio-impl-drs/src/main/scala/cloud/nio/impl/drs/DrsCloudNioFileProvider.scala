package cloud.nio.impl.drs

import java.nio.channels.{ReadableByteChannel, WritableByteChannel}

import cloud.nio.spi.{CloudNioFileList, CloudNioFileProvider, CloudNioRegularFileAttributes}
import com.typesafe.config.Config
import org.apache.http.HttpStatus
import org.apache.http.client.methods.CloseableHttpResponse
import org.apache.http.impl.client.{CloseableHttpClient, HttpClientBuilder}

import scala.util.Try


class DrsCloudNioFileProvider(config: Config, scheme: String, drsPathResolver: DrsPathResolver) extends CloudNioFileProvider {

  private def getDrsPath(cloudHost: String, cloudPath: String): String = s"$scheme://$cloudHost/$cloudPath"

  private def checkIfPathExistsThroughMartha(drsPath: String): Boolean = {
    val httpClient: CloseableHttpClient = HttpClientBuilder.create().build()

    try {
      val marthaResponse: CloseableHttpResponse = drsPathResolver.makeHttpRequestToMartha(drsPath, httpClient)

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


  override def existsPath(cloudHost: String, cloudPath: String): Boolean =
    checkIfPathExistsThroughMartha(getDrsPath(cloudHost, cloudPath))

  //This method is used to check if the given path is a directory or not. Since DRS currently doesn't
  //support resolving url to directories, it is false
  override def existsPaths(cloudHost: String, cloudPathPrefix: String): Boolean = false

  override def listObjects(cloudHost: String, cloudPathPrefix: String, markerOption: Option[String]): CloudNioFileList =
    throw new UnsupportedOperationException("DRS currently doesn't support listObjects.")

  override def copy(sourceCloudHost: String, sourceCloudPath: String, targetCloudHost: String, targetCloudPath: String): Unit =
    throw new UnsupportedOperationException("DRS currently doesn't support copy.")

  override def deleteIfExists(cloudHost: String, cloudPath: String): Boolean =
    throw new UnsupportedOperationException("DRS currently doesn't support delete.")

  override def read(cloudHost: String, cloudPath: String, offset: Long): ReadableByteChannel =
    throw new UnsupportedOperationException("DRS currently doesn't support read.")

  override def write(cloudHost: String, cloudPath: String): WritableByteChannel =
    throw new UnsupportedOperationException("DRS currently doesn't support write.")

  override def fileAttributes(cloudHost: String, cloudPath: String): Option[CloudNioRegularFileAttributes] =
    Option(new DrsCloudNioRegularFileAttributes(config, getDrsPath(cloudHost,cloudPath), drsPathResolver))
}
