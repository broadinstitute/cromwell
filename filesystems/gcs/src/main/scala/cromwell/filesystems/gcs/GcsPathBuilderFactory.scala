package cromwell.filesystems.gcs

import akka.actor.ActorSystem
import com.google.api.client.googleapis.media.MediaHttpUploader
import com.google.api.gax.retrying.RetrySettings
import com.google.auth.Credentials
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import cromwell.filesystems.gcs.auth.GoogleAuthMode
import net.ceedubs.ficus.Ficus._

import scala.concurrent.ExecutionContext

object GcsPathBuilderFactory {

  private[this] lazy val UploadBufferBytes = {
    ConfigFactory.load().as[Option[Int]]("google.upload-buffer-bytes").getOrElse(MediaHttpUploader.MINIMUM_CHUNK_SIZE)
  }

  val DefaultCloudStorageConfiguration = {
    CloudStorageConfiguration.builder()
      .blockSize(UploadBufferBytes)
      .permitEmptyPathComponents(true)
      .stripPrefixSlash(true)
      .usePseudoDirectories(true)
      .build()
  }
}

case class GcsPathBuilderFactory(authMode: GoogleAuthMode,
                                 applicationName: String,
                                 retrySettings: Option[RetrySettings] = None,
                                 cloudStorageConfiguration: CloudStorageConfiguration = GcsPathBuilderFactory.DefaultCloudStorageConfiguration)

  extends PathBuilderFactory {

  def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext) = {
    GcsPathBuilder.fromAuthMode(authMode, applicationName, retrySettings, cloudStorageConfiguration, options)
  }

  /**
    * Ignores the authMode and creates a GcsPathBuilder using the passed credentials directly.
    * Can be used when the Credentials are already available.
    */
  def fromCredentials(options: WorkflowOptions, credentials: Credentials) = {
    GcsPathBuilder.fromCredentials(credentials, applicationName, retrySettings, cloudStorageConfiguration, options)
  }
}
