package cromwell.filesystems.gcs

import java.net.URI

import akka.actor.ActorSystem
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.gax.retrying.RetrySettings
import com.google.auth.Credentials
import com.google.cloud.http.HttpTransportOptions
import com.google.cloud.storage.contrib.nio.{CloudStorageConfiguration, CloudStorageFileSystem, CloudStorageFileSystemProvider, CloudStoragePath}
import com.google.cloud.storage.{BlobId, StorageOptions}
import com.google.common.base.Preconditions._
import com.google.common.net.UrlEscapers
import cromwell.core.WorkflowOptions
import cromwell.core.path.{NioPath, Path, PathBuilder}
import cromwell.filesystems.gcs.GcsPathBuilder._
import cromwell.filesystems.gcs.auth.GoogleAuthMode

import scala.concurrent.{ExecutionContext, Future}
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.Try

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

  sealed trait GcsPathValidation
  case object ValidFullGcsPath extends GcsPathValidation
  case object PossiblyValidRelativeGcsPath extends GcsPathValidation
  sealed trait InvalidGcsPath extends GcsPathValidation {
    def pathString: String
    def errorMessage: String
  }
  final case class InvalidFullGcsPath(pathString: String) extends InvalidGcsPath {
    override def errorMessage = {
      val prefix = s"""
      |The bucket name in GCS path '$pathString' is not compatible with URI host name standards.
      |URI host name compatibility is a requirement for Cromwell's GCS filesystem support.
      |Google also generally advises against the use of underscores in GCS bucket names, as well as against
      |the use of periods or dashes in certain patterns as described here:
      |https://cloud.google.com/storage/docs/naming.
      """.stripMargin.replaceAll("\n", " ").trim
      val underscoreWarning = if (pathString.contains("_")) s"In particular, the bucket name in '$pathString' may contain an underscore which is not a valid character in a URI host." else ""
      List(prefix, underscoreWarning).mkString(" ")
    }
  }
  final case class UnparseableGcsPath(pathString: String, throwable: Throwable) extends InvalidGcsPath {
    override def errorMessage: String =
      List(s"The specified GCS path '$pathString' does not parse as a URI.", throwable.getMessage).mkString("\n")
  }

  def validateGcsPath(string: String): GcsPathValidation = {
    Try {
      val uri = getUri(string)
      if (uri.getScheme == null) PossiblyValidRelativeGcsPath
      else if (uri.getScheme == "gs") {
        if (uri.getHost == null) InvalidFullGcsPath(string) else ValidFullGcsPath
      } else InvalidFullGcsPath(string)
    } recover { case t => UnparseableGcsPath(string, t) } get
  }

  def isGcsPath(nioPath: NioPath): Boolean = {
    nioPath.getFileSystem.provider().getScheme == CloudStorageFileSystem.URI_SCHEME
  }

  def getUri(string: String) = URI.create(UrlEscapers.urlFragmentEscaper().escape(string))
  
  def fromAuthMode(authMode: GoogleAuthMode,
                   applicationName: String,
                   retrySettings: Option[RetrySettings],
                   cloudStorageConfiguration: CloudStorageConfiguration,
                   options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[GcsPathBuilder] = {
    authMode.credential(options) map { credentials =>
      fromCredentials(credentials,
        applicationName,
        retrySettings,
        cloudStorageConfiguration,
        options
      )
    }
  }
  
  def fromCredentials(credentials: Credentials,
                      applicationName: String,
                      retrySettings: Option[RetrySettings],
                      cloudStorageConfiguration: CloudStorageConfiguration,
                      options: WorkflowOptions): GcsPathBuilder = {
    val transportOptions = HttpTransportOptions.newBuilder()
      .setReadTimeout(3.minutes.toMillis.toInt)
      .build()

    val storageOptionsBuilder = StorageOptions.newBuilder()
      .setTransportOptions(transportOptions)
      .setCredentials(credentials)

    retrySettings foreach storageOptionsBuilder.setRetrySettings

    // Grab the google project from Workflow Options if specified and set
    // that to be the project used by the StorageOptions Builder
    options.get("google_project") map storageOptionsBuilder.setProjectId


    val storageOptions = storageOptionsBuilder.build()

    // Create a com.google.api.services.storage.Storage
    // This is the underlying api used by com.google.cloud.storage
    // By bypassing com.google.cloud.storage, we can create low level requests that can be batched
    val apiStorage: com.google.api.services.storage.Storage = {
      new com.google.api.services.storage.Storage
      .Builder(HttpTransport, JsonFactory, GoogleConfiguration.withCustomTimeouts(transportOptions.getHttpRequestInitializer(storageOptions)))
        .setApplicationName(applicationName)
        .build()
    }

    // Create a com.google.cloud.storage.Storage
    // This is the "relatively" high level API, and recommended one. The nio implementation sits on top of this.
    val cloudStorage: com.google.cloud.storage.Storage = storageOptions.getService

    // The CloudStorageFileSystemProvider constructor is not public. Currently the only way to obtain one is through a CloudStorageFileSystem
    // Moreover at this point we can use the same provider for all operations as we have usable credentials
    // In order to avoid recreating a provider with every getPath call, create a dummy FileSystem just to get its provider
    val provider: CloudStorageFileSystemProvider = CloudStorageFileSystem.forBucket("dummy", cloudStorageConfiguration, storageOptions).provider()
    
    new GcsPathBuilder(apiStorage, cloudStorage, provider, storageOptions.getProjectId)
  }
}

class GcsPathBuilder(val apiStorage: com.google.api.services.storage.Storage,
                     val cloudStorage: com.google.cloud.storage.Storage,
                     provider: CloudStorageFileSystemProvider,
                     val projectId: String) extends PathBuilder {
  def build(string: String): Try[GcsPath] = {
    Try {
      val uri = getUri(string)
      GcsPathBuilder.checkValid(uri)
      GcsPath(provider.getPath(uri), apiStorage, cloudStorage)
    }
  }

  override def name: String = "Gcs"
}

case class GcsPath private[gcs](nioPath: NioPath,
                                apiStorage: com.google.api.services.storage.Storage,
                                cloudStorage: com.google.cloud.storage.Storage
                               ) extends Path {
  lazy val blob = BlobId.of(cloudStoragePath.bucket, cloudStoragePath.toRealPath().toString)

  override protected def newPath(nioPath: NioPath): GcsPath = GcsPath(nioPath, apiStorage, cloudStorage)

  override def pathAsString: String = java.net.URLDecoder.decode(nioPath.toUri.toString, "UTF-8")

  override def pathWithoutScheme: String = {
    val gcsPath = cloudStoragePath
    gcsPath.bucket + gcsPath.toAbsolutePath.toString
  }

  def cloudStoragePath: CloudStoragePath = nioPath match {
    case gcsPath: CloudStoragePath => gcsPath
    case _ => throw new RuntimeException(s"Internal path was not a cloud storage path: $nioPath")
  }
}
