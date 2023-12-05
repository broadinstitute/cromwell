package cromwell.cloudsupport.gcp.gcs

import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.googleapis.media.MediaHttpUploader
import com.google.api.client.json.gson.GsonFactory
import com.google.api.gax.retrying.RetrySettings
import com.google.cloud.storage.StorageOptions
import com.google.api.services.storage.Storage
import com.google.auth.Credentials
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration
import com.typesafe.config.ConfigFactory
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.http.GoogleHttpTransportOptions.TransportOptions
import net.ceedubs.ficus.Ficus._

object GcsStorage {
  val JsonFactory = GsonFactory.getDefaultInstance
  val HttpTransport = GoogleNetHttpTransport.newTrustedTransport

  val DefaultCloudStorageConfiguration = {
    val UploadBufferBytes =
      ConfigFactory.load().as[Option[Int]]("google.upload-buffer-bytes").getOrElse(MediaHttpUploader.MINIMUM_CHUNK_SIZE)

    CloudStorageConfiguration
      .builder()
      .blockSize(UploadBufferBytes)
      .permitEmptyPathComponents(true)
      .stripPrefixSlash(true)
      .usePseudoDirectories(true)
      .build()
  }

  def gcsStorage(applicationName: String, storageOptions: StorageOptions): Storage =
    new Storage.Builder(
      HttpTransport,
      JsonFactory,
      GoogleConfiguration.withCustomTimeouts(TransportOptions.getHttpRequestInitializer(storageOptions))
    )
      .setApplicationName(applicationName)
      .build()

  def gcsStorage(applicationName: String, credentials: Credentials, retrySettings: RetrySettings): Storage =
    gcsStorage(applicationName, gcsStorageOptions(credentials, retrySettings))

  def gcsStorageOptions(credentials: Credentials,
                        retrySettings: RetrySettings,
                        project: Option[String] = None
  ): StorageOptions = {
    val storageOptionsBuilder = StorageOptions
      .newBuilder()
      .setTransportOptions(TransportOptions)
      .setCredentials(credentials)

    storageOptionsBuilder.setRetrySettings(retrySettings)
    project foreach storageOptionsBuilder.setProjectId

    storageOptionsBuilder.build()
  }
}
