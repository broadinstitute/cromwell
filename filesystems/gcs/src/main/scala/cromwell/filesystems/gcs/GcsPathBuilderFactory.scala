package cromwell.filesystems.gcs

import akka.actor.ActorSystem
import com.google.api.client.googleapis.media.MediaHttpUploader
import com.google.api.gax.retrying.RetrySettings
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import cromwell.filesystems.gcs.auth.GoogleAuthMode
import net.ceedubs.ficus.Ficus._

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

  def withOptions(options: WorkflowOptions)(implicit actorSystem: ActorSystem) = new GcsPathBuilder(authMode, applicationName, retrySettings, cloudStorageConfiguration, options)
}
