package cromwell.filesystems.gcs

import java.net.URI
import java.nio.file.Path

import akka.actor.ActorSystem
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.cloud.RetryParams
import com.google.cloud.storage.StorageOptions
import com.google.cloud.storage.contrib.nio.{CloudStorageConfiguration, CloudStorageFileSystem, CloudStoragePath}
import com.google.common.base.Preconditions._
import cromwell.core.WorkflowOptions
import cromwell.core.path.{CustomRetryParams, PathBuilder, RetryableFileSystemProviderProxy}
import cromwell.filesystems.gcs.auth.GoogleAuthMode

import scala.util.Try

object GcsPathBuilder {

  val JsonFactory = JacksonFactory.getDefaultInstance
  val HttpTransport = GoogleNetHttpTransport.newTrustedTransport

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
  import GcsPathBuilder._

  authMode.validate(options)

  protected val storageOptions = StorageOptions.builder()
    .authCredentials(authMode.authCredentials(options))
    .retryParams(retryParams)
    .build()

  // The CloudStorageFileSystemProvider constructor is not public. Currently the only way to obtain one is through a CloudStorageFileSystem
  // Moreover at this point we can use the same provider for all operations as we have usable credentials
  // In order to avoid recreating a provider with every getPath call, create a dummy FileSystem just to get its provider
  protected val provider = CloudStorageFileSystem.forBucket("dummy", cloudStorageConfiguration, storageOptions).provider()

  /*
   * The StorageService already contains a StorageRpc object that contains a com.google.api.services.storage.Storage object
   * However it is not accessible from StorageService.
   * com.google.cloud.storage.Storage has some batching capabilities but not for copying.
   * In order to support batch copy, we need a com.google.api.services.storage.Storage.
   */
  private lazy val apiStorage = new com.google.api.services.storage.Storage(HttpTransport, JsonFactory, authMode.credential(options))

  def getHash(path: CloudStoragePath) = {
    storageOptions.service().get(path.bucket(), path.toRealPath().toString).crc32c()
  }

  def batchCopy(pairs: Seq[(Path, Path)]) = {
    // TODO: implement with apiStorage
  }

  def build(string: String): Try[Path] = {
    Try {
      val uri = URI.create(string)
      GcsPathBuilder.checkValid(uri)
      provider.getPath(uri)
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

  override protected val provider = new RetryableFileSystemProviderProxy(super.provider, customRetryParams)

  override def getHash(path: CloudStoragePath) = provider.withRetry(() => super.getHash(path))
}
