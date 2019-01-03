package cloud.nio.impl.drs

import java.io.ByteArrayInputStream
import java.nio.channels.{ReadableByteChannel, WritableByteChannel}

import cloud.nio.spi.{CloudNioFileList, CloudNioFileProvider, CloudNioRegularFileAttributes}
import com.google.auth.oauth2.GoogleCredentials
import com.google.cloud.ReadChannel
import com.google.cloud.storage.StorageOptions
import com.typesafe.config.Config
import org.apache.http.HttpStatus
import org.apache.http.client.methods.CloseableHttpResponse
import org.apache.http.impl.client.{CloseableHttpClient, HttpClientBuilder}

import scala.util.Try


class DrsCloudNioFileProvider(config: Config, scheme: String, drsPathResolver: DrsPathResolver) extends CloudNioFileProvider {

  private lazy val httpClientBuilder = HttpClientBuilder.create()

  private val GcsScheme = "gs"


  private def getDrsPath(cloudHost: String, cloudPath: String): String = s"$scheme://$cloudHost/$cloudPath"


  private def checkIfPathExistsThroughMartha(drsPath: String): Boolean = {
    val httpClient: CloseableHttpClient = httpClientBuilder.build()

    try {
      val marthaResponse: CloseableHttpResponse = drsPathResolver.makeHttpRequestToMartha(drsPath, httpClient, needServiceAccount = false)

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


  private def extractUrlForScheme(drsPath: String, urlArray: Array[Url], scheme: String): String = {
    val schemeUrlOption = urlArray.find(_.url.startsWith(scheme))

    schemeUrlOption match {
      case Some(schemeUrl) => schemeUrl.url
      case None => throw UrlNotFoundException(drsPath, scheme)
    }
  }


  private def gcsInputStream(bucket: String, fileToBeLocalized: String, serviceAccount: String): ReadChannel = {
    val credentials = GoogleCredentials.fromStream(new ByteArrayInputStream(serviceAccount.getBytes()))
    val storage = StorageOptions.newBuilder().setCredentials(credentials).build().getService
    val blob = storage.get(bucket, fileToBeLocalized)
    blob.reader()
  }


  private def inputStream(url: String, serviceAccount: String, scheme: String): ReadChannel = {
    scheme match {
      case GcsScheme => {
        val Array(bucket, fileToBeLocalized) = url.replace(s"$scheme://", "").split("/", 2)
        gcsInputStream(bucket, fileToBeLocalized, serviceAccount)
      }
      case otherScheme => throw new UnsupportedOperationException(s"DRS currently doesn't support reading files for $otherScheme.")
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
    val marthaResponse = drsPathResolver.resolveDrsThroughMartha(getDrsPath(cloudHost,cloudPath), needServiceAccount = true)

    val serviceAccount = marthaResponse.googleServiceAccount match {
      case Some(googleSA) => googleSA.data.toString
      case None => throw GoogleSANotFoundException(getDrsPath(cloudHost, cloudPath))
    }

    //Currently, Martha only supports resolving DRS paths to GCS paths
    val url = extractUrlForScheme(getDrsPath(cloudHost, cloudPath), marthaResponse.dos.data_object.urls, GcsScheme)
    inputStream(url, serviceAccount, GcsScheme)
  }


  override def write(cloudHost: String, cloudPath: String): WritableByteChannel =
    throw new UnsupportedOperationException("DRS currently doesn't support write.")


  override def fileAttributes(cloudHost: String, cloudPath: String): Option[CloudNioRegularFileAttributes] =
    Option(new DrsCloudNioRegularFileAttributes(config, getDrsPath(cloudHost,cloudPath), drsPathResolver))
}


case class GoogleSANotFoundException(drsPath: String) extends Exception(s"Error finding Google Service Account associated with DRS path $drsPath through Martha.")

case class UrlNotFoundException(drsPath: String, scheme: String) extends Exception(s"DRS was not able to find a $scheme url associated with $drsPath.")

