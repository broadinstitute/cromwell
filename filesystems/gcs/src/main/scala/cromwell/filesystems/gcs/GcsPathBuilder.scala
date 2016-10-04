package cromwell.filesystems.gcs

import java.net.URI
import java.nio.file.Path

import akka.actor.ActorSystem
import com.google.cloud.RetryParams
import com.google.cloud.storage.StorageOptions
import com.google.cloud.storage.contrib.nio.{CloudStorageConfiguration, CloudStorageFileSystem, CloudStoragePath}
import com.google.common.base.Preconditions._
import cromwell.core.WorkflowOptions
import cromwell.core.path.{CustomRetryParams, PathBuilder, RetryableFileSystemProviderProxy}
import cromwell.filesystems.gcs.auth.GoogleAuthMode

import scala.util.Try

object GcsPathBuilder {

  def checkValid(uri: URI) = {
    checkNotNull(uri.getScheme, s"%s does not have a gcs scheme", uri)
    checkArgument(
      uri.getScheme.equalsIgnoreCase(CloudStorageFileSystem.URI_SCHEME),
      "Cloud Storage URIs must have '%s' scheme: %s",
      CloudStorageFileSystem.URI_SCHEME,
      uri
    )
    checkNotNull(uri.getHost, s"%s does not have a host", uri)
  }

  def isValidGcsUrl(str: String): Boolean = {
    Try(checkValid(URI.create(str))).isSuccess
  }
}

class GcsPathBuilder(authMode: GoogleAuthMode,
                     retryParams: RetryParams,
                     cloudStorageConfiguration: CloudStorageConfiguration,
                     options: WorkflowOptions) extends PathBuilder {

  authMode.validate(options)

  protected val storageOptions = StorageOptions.builder()
    .authCredentials(authMode.authCredentials(options))
    .retryParams(retryParams)
    .build()

  def getHash(path: CloudStoragePath) = {
    storageOptions.service().get(path.bucket(), path.toRealPath().toString).crc32c()
  }

  def build(string: String): Try[Path] = {
    Try {
      val uri = URI.create(string)
      GcsPathBuilder.checkValid(uri)
      val gcsFileSystem = CloudStorageFileSystem.forBucket(uri.getHost, cloudStorageConfiguration, storageOptions)
      gcsFileSystem.getPath(uri.getPath)
    }
  }

  override def name: String = "Gcs"
}

class RetryableGcsPathBuilder(authMode: GoogleAuthMode,
                              googleRetryParams: RetryParams,
                              customRetryParams: CustomRetryParams,
                              cloudStorageConfiguration: CloudStorageConfiguration,
                              options: WorkflowOptions)(implicit actorSystem: ActorSystem)
  extends GcsPathBuilder(authMode, googleRetryParams, cloudStorageConfiguration, options) {

  override def build(string: String): Try[CloudStoragePath] = Try {
    val uri = URI.create(string)
    GcsPathBuilder.checkValid(uri)
    val gcsFileSystem = CloudStorageFileSystem.forBucket(uri.getHost, cloudStorageConfiguration, storageOptions)
    new RetryableFileSystemProviderProxy(gcsFileSystem.provider(), customRetryParams).getPath(uri).asInstanceOf[CloudStoragePath]
  }
}
