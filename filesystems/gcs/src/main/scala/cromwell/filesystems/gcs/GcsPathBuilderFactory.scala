package cromwell.filesystems.gcs

import akka.actor.ActorSystem
import com.google.api.client.googleapis.media.MediaHttpUploader
import com.google.cloud.RetryParams
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import cromwell.filesystems.gcs.auth.GoogleAuthMode
import lenthall.config.ScalaConfig._

object GcsPathBuilderFactory {

  private[this] lazy val UploadBufferBytes = {
    ConfigFactory.load().getBytesOr("google.upload-buffer-bytes", MediaHttpUploader.MINIMUM_CHUNK_SIZE).toInt
  }

  val DefaultRetryParams = RetryParams.defaultInstance()
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
                                 retryParams: RetryParams = GcsPathBuilderFactory.DefaultRetryParams,
                                 cloudStorageConfiguration: CloudStorageConfiguration = GcsPathBuilderFactory.DefaultCloudStorageConfiguration)

  extends PathBuilderFactory {

  def withOptions(options: WorkflowOptions) = new GcsPathBuilder(authMode, retryParams, cloudStorageConfiguration, options)
}

case class RetryableGcsPathBuilderFactory(authMode: GoogleAuthMode,
                                 retryParams: RetryParams = GcsPathBuilderFactory.DefaultRetryParams,
                                 cloudStorageConfiguration: CloudStorageConfiguration = GcsPathBuilderFactory.DefaultCloudStorageConfiguration)
                                         (implicit actorSystem: ActorSystem)

  extends PathBuilderFactory {

  def withOptions(options: WorkflowOptions) = new RetryableGcsPathBuilder(authMode, retryParams, cloudStorageConfiguration, options)
}
