package cromwell.filesystems.gcs

import java.net.URI
import java.nio.file.Path

import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.cloud.RetryParams
import com.google.cloud.storage.StorageOptions
import com.google.cloud.storage.contrib.nio.{CloudStorageConfiguration, CloudStorageFileSystem, CloudStoragePath}
import com.google.common.base.Preconditions._
import com.google.common.net.UrlEscapers
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilder
import cromwell.filesystems.gcs.GcsPathBuilder._
import cromwell.filesystems.gcs.auth.GoogleAuthMode

import scala.util.{Failure, Try}

object GcsPathBuilder {

  val JsonFactory = JacksonFactory.getDefaultInstance
  val HttpTransport = GoogleNetHttpTransport.newTrustedTransport

  def checkValid(uri: URI) = {
    checkNotNull(uri.getScheme, s"%s does not have a gcs scheme", uri)
    checkArgument(
      uri.getScheme.equalsIgnoreCase(CloudStorageFileSystem.URI_SCHEME),
      "Cloud Storage URIs must have '%s' scheme: %s",
      CloudStorageFileSystem.URI_SCHEME: Any,
      uri: Any
    )
    checkNotNull(uri.getHost, s"%s does not have a host", uri)
  }

  def isValidGcsUrl(str: String): Boolean = {
    Try(checkValid(getUri(str))).isSuccess
  }

  def isGcsPath(path: Path): Boolean = {
    path.getFileSystem.provider().getScheme == CloudStorageFileSystem.URI_SCHEME
  }

  def getUri(string: String) = URI.create(UrlEscapers.urlFragmentEscaper().escape(string))
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

  // The CloudStorageFileSystemProvider constructor is not public. Currently the only way to obtain one is through a CloudStorageFileSystem
  // Moreover at this point we can use the same provider for all operations as we have usable credentials
  // In order to avoid recreating a provider with every getPath call, create a dummy FileSystem just to get its provider
  private val provider = CloudStorageFileSystem.forBucket("dummy", cloudStorageConfiguration, storageOptions).provider()

  /*
   * The StorageService already contains a StorageRpc object that contains a com.google.api.services.storage.Storage object
   * However it is not accessible from StorageService.
   * com.google.cloud.storage.Storage has some batching capabilities but not for copying.
   * In order to support batch copy, we need a com.google.api.services.storage.Storage.
   */
  def getHash(path: Path): Try[String] = {
    path match {
      case gcsPath: CloudStoragePath => Try(storageOptions.service().get(gcsPath.bucket(), gcsPath.toRealPath().toString).crc32c())
      case other => Failure(new IllegalArgumentException(s"$other is not a CloudStoragePath"))
    }
  }

  def build(string: String): Try[Path] = {
    Try {
      val uri = getUri(string)
      GcsPathBuilder.checkValid(uri)
      provider.getPath(uri)
    }
  }

  override def name: String = "Gcs"
}
