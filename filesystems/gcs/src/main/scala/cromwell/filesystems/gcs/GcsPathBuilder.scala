package cromwell.filesystems.gcs

import java.io._
import java.net.URI
import java.nio.channels.Channels
import java.nio.file.NoSuchFileException

import akka.actor.ActorSystem
import akka.http.scaladsl.model.{ContentTypes, StatusCodes}
import better.files.File.OpenOptions
import cats.effect.IO
import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.google.api.gax.retrying.RetrySettings
import com.google.auth.Credentials
import com.google.cloud.storage.Storage.{BlobField, BlobGetOption, BlobSourceOption, BlobTargetOption}
import com.google.cloud.storage.contrib.nio.{CloudStorageConfiguration, CloudStorageFileSystem, CloudStoragePath}
import com.google.cloud.storage.{BlobId, BlobInfo, StorageOptions}
import com.google.common.cache.CacheBuilder
import com.google.common.net.UrlEscapers
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.cloudsupport.gcp.gcs.GcsStorage
import cromwell.core.WorkflowOptions
import cromwell.core.path.{NioPath, Path, PathBuilder}
import cromwell.filesystems.gcs.GcsPath.BucketInformationGetter
import cromwell.filesystems.gcs.GcsPathBuilder._
import cromwell.filesystems.gcs.GoogleUtil._
import cromwell.filesystems.gcs.cache.GcsBucketCache.GcsBucketInformation
import cromwell.filesystems.gcs.cache.GcsBucketInformationPolicies.GcsBucketInformationPolicy
import cromwell.filesystems.gcs.cache.GcsFileSystemCache
import mouse.all._

import scala.concurrent.{ExecutionContext, Future}
import scala.io.Codec
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
object GcsPathBuilder {
  implicit class EnhancedCromwellPath(val path: Path) extends AnyVal {
    def requesterPaysProject: Option[String] = path match {
      case gcs: GcsPath => gcs.requesterPays.option(gcs.projectId)
      case _ => None
    }

    def requesterPaysGSUtilFlagList = requesterPaysProject.map(project => List("-u", project)).getOrElse(List.empty)
    def requesterPaysGSUtilFlag = requesterPaysGSUtilFlagList.mkString(" ")
  }

  /*
    * Provides some level of validation of GCS bucket names
    * This is meant to alert the user early if they mistyped a gcs path in their workflow / inputs and not to validate
    * exact bucket syntax, which is done by GCS.
    * See https://cloud.google.com/storage/docs/naming for full spec
  */
  val GcsBucketPattern =
    """
      (?x)                                      # Turn on comments and whitespace insensitivity
      ^gs://
      (                                         # Begin capturing group for gcs bucket name
        [a-z0-9][a-z0-9-_\\.]+[a-z0-9]          # Regex for bucket name - soft validation, see comment above
      )                                         # End capturing group for gcs bucket name
      (?:
        /.*                                     # No validation here
      )?
    """.trim.r

  sealed trait GcsPathValidation
  case class ValidFullGcsPath(bucket: String, path: String) extends GcsPathValidation
  case object PossiblyValidRelativeGcsPath extends GcsPathValidation
  sealed trait InvalidGcsPath extends GcsPathValidation {
    def pathString: String
    def errorMessage: String
  }
  final case class InvalidScheme(pathString: String) extends InvalidGcsPath {
    override def errorMessage = s"Cloud Storage URIs must have 'gs' scheme: $pathString"
  }
  final case class InvalidFullGcsPath(pathString: String) extends InvalidGcsPath {
    override def errorMessage = {
      s"""
         |The path '$pathString' does not seem to be a valid GCS path.
         |Please check that it starts with gs:// and that the bucket and object follow GCS naming guidelines at
         |https://cloud.google.com/storage/docs/naming.
      """.stripMargin.replace("\n", " ").trim
    }
  }
  final case class UnparseableGcsPath(pathString: String, throwable: Throwable) extends InvalidGcsPath {
    override def errorMessage: String =
      List(s"The specified GCS path '$pathString' does not parse as a URI.", throwable.getMessage).mkString("\n")
  }

  /**
    * Tries to extract a bucket name out of the provided string using rules less strict that URI hostname,
    * as GCS allows albeit discourages.
    */
  private def softBucketParsing(string: String): Option[String] = string match {
    case GcsBucketPattern(bucket) => Option(bucket)
    case _ => None
  }

  def validateGcsPath(string: String): GcsPathValidation = {
    Try {
      val uri = URI.create(UrlEscapers.urlFragmentEscaper().escape(string))
      if (uri.getScheme == null) PossiblyValidRelativeGcsPath
      else if (uri.getScheme.equalsIgnoreCase(CloudStorageFileSystem.URI_SCHEME)) {
        if (uri.getHost == null) {
          softBucketParsing(string) map { ValidFullGcsPath(_, uri.getPath) } getOrElse InvalidFullGcsPath(string)
        } else ValidFullGcsPath(uri.getHost, uri.getPath)
      } else InvalidScheme(string)
    } recover { case t => UnparseableGcsPath(string, t) } get
  }

  def isGcsPath(nioPath: NioPath): Boolean = {
    nioPath.getFileSystem.provider().getScheme.equalsIgnoreCase(CloudStorageFileSystem.URI_SCHEME)
  }

  def fromAuthMode(authMode: GoogleAuthMode,
                   applicationName: String,
                   retrySettings: RetrySettings,
                   cloudStorageConfiguration: CloudStorageConfiguration,
                   options: WorkflowOptions,
                   defaultProject: Option[String],
                   bucketInformationPolicy: GcsBucketInformationPolicy)(implicit as: ActorSystem, ec: ExecutionContext): Future[GcsPathBuilder] = {
    authMode.retryCredential(options) map { credentials =>
      fromCredentials(credentials,
        applicationName,
        retrySettings,
        cloudStorageConfiguration,
        options,
        defaultProject,
        bucketInformationPolicy
      )
    }
  }

  def fromCredentials(credentials: Credentials,
                      applicationName: String,
                      retrySettings: RetrySettings,
                      cloudStorageConfiguration: CloudStorageConfiguration,
                      options: WorkflowOptions,
                      defaultProject: Option[String],
                      bucketInformationPolicy: GcsBucketInformationPolicy)(implicit ec: ExecutionContext): GcsPathBuilder = {
    // Grab the google project from Workflow Options if specified and set
    // that to be the project used by the StorageOptions Builder. If it's not
    // specified use the default project mentioned in config file
    val project: Option[String] =  options.get("google_project").toOption match {
      case Some(googleProject) => Option(googleProject)
      case None => defaultProject
    }

    val storageOptions = GcsStorage.gcsStorageOptions(credentials, retrySettings, project)

    // Create a com.google.api.services.storage.Storage
    // This is the underlying api used by com.google.cloud.storage
    // By bypassing com.google.cloud.storage, we can create low level requests that can be batched
    val apiStorage = GcsStorage.gcsStorage(applicationName, storageOptions)

    new GcsPathBuilder(apiStorage, cloudStorageConfiguration, storageOptions, bucketInformationPolicy)
  }
}

class GcsPathBuilder(apiStorage: com.google.api.services.storage.Storage,
                     cloudStorageConfiguration: CloudStorageConfiguration,
                     storageOptions: StorageOptions,
                     bucketInformationPolicy: GcsBucketInformationPolicy)(implicit ec: ExecutionContext) extends PathBuilder {
  private [gcs] val projectId = storageOptions.getProjectId
  private [gcs] val bucketStore = bucketInformationPolicy.toStore(cloudStorage, projectId)
  // We can cache filesystems here since they only depend on the bucket and the credentials (which are unique per instance of path builder)
  private val fileSystemCache = new GcsFileSystemCache(cloudStorage, CacheBuilder.newBuilder().build[String, CloudStorageFileSystem](), cloudStorageConfiguration, storageOptions)

  private lazy val cloudStorage = storageOptions.getService

  /**
    * Tries to create a new GcsPath from a String representing an absolute gcs path: gs://<bucket>[/<path>].
    *
    * Note that this creates a new CloudStorageFileSystemProvider for every Path created, hence making it unsuitable for
    * file copying using the nio copy method. Cromwell currently uses a lower level API which works around this problem.
    *
    * If you plan on using the nio copy method make sure to take this into consideration.
    *
    * Also see https://github.com/GoogleCloudPlatform/google-cloud-java/issues/1343
    */
  def build(string: String): Try[GcsPath] = {
    validateGcsPath(string) match {
      case ValidFullGcsPath(bucket, path) =>
        val pathAndRpIO = for {
          fileSystem <- fileSystemCache.getCachedValue(bucket)
          cloudStoragePath = fileSystem.getPath(path)
        } yield cloudStoragePath
        
        Try {
          val cloudStoragePath = pathAndRpIO.unsafeRunSync()
          GcsPath(cloudStoragePath, apiStorage, cloudStorage, projectId, bucketInformation)
        }
      case PossiblyValidRelativeGcsPath => Failure(new IllegalArgumentException(s"$string does not have a gcs scheme"))
      case invalid: InvalidGcsPath => Failure(new IllegalArgumentException(invalid.errorMessage))
    }
  }
  
  // Provide a way to the Path to get most up to date information about the status of requester pays for its bucket
  private val bucketInformation: BucketInformationGetter = bucket => for {
    // Use the provided EC to get bucket information
    _ <- IO.shift(ec)
    bucketInfo <- bucketStore.information(bucket)
  } yield bucketInfo

  override def name: String = "Google Cloud Storage"
}

object GcsPath {
  type BucketInformationGetter = String => IO[GcsBucketInformation]
}

case class GcsPath private[gcs](nioPath: NioPath,
                                apiStorage: com.google.api.services.storage.Storage,
                                cloudStorage: com.google.cloud.storage.Storage,
                                projectId: String,
                                bucketInformationGetter: BucketInformationGetter) extends Path {
  lazy val blob = BlobId.of(cloudStoragePath.bucket, cloudStoragePath.toRealPath().toString)

  def requesterPays = bucketInformationGetter(blob.getBucket).unsafeRunSync().requesterPays
  /**
    * Will be null if requesterPays is false 
    */
  def requesterPaysProject = requesterPays.option(projectId).orNull
  // Utility flags to pass to gcs requests
  private def userProjectBlobTarget: List[BlobTargetOption] = requesterPays.option(BlobTargetOption.userProject(projectId)).toList
  private def userProjectBlobSource: List[BlobSourceOption] = requesterPays.option(BlobSourceOption.userProject(projectId)).toList
  private def userProjectBlobGet: List[BlobGetOption] = requesterPays.option(BlobGetOption.userProject(projectId)).toList

  override protected def newPath(nioPath: NioPath): GcsPath = GcsPath(nioPath, apiStorage, cloudStorage, projectId, bucketInformationGetter)

  override def pathAsString: String = {
    val host = cloudStoragePath.bucket().stripSuffix("/")
    val path = cloudStoragePath.toString.stripPrefix("/")
    s"${CloudStorageFileSystem.URI_SCHEME}://$host/$path"
  }

  override def writeContent(content: String)(openOptions: OpenOptions, codec: Codec) = {
    cloudStorage.create(
      BlobInfo.newBuilder(blob)
        .setContentType(ContentTypes.`text/plain(UTF-8)`.value)
        .build(),
      content.getBytes(codec.charSet),
      userProjectBlobTarget: _*
    )
    this
  }

  override def mediaInputStream: InputStream = {
    Try{
      Channels.newInputStream(cloudStorage.reader(blob, userProjectBlobSource: _*))
    } match {
      case Success(inputStream) => inputStream
      case Failure(e: GoogleJsonResponseException) if e.getStatusCode == StatusCodes.NotFound.intValue =>
        throw new NoSuchFileException(pathAsString)
      case Failure(e) => e.getMessage
        throw new IOException(s"Failed to open an input stream for $pathAsString: ${e.getMessage}", e)
    }
  }

  override def hash = {
    cloudStorage.get(blob, userProjectBlobGet :+ BlobGetOption.fields(BlobField.CRC32C): _*).getCrc32c
  }

  override def pathWithoutScheme: String = {
    val gcsPath = cloudStoragePath
    gcsPath.bucket + gcsPath.toAbsolutePath.toString
  }

  def cloudStoragePath: CloudStoragePath = nioPath match {
    case gcsPath: CloudStoragePath => gcsPath
    case _ => throw new RuntimeException(s"Internal path was not a cloud storage path: $nioPath")
  }
}
