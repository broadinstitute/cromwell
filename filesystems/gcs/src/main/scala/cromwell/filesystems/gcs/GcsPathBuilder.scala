package cromwell.filesystems.gcs

import akka.actor.ActorSystem
import akka.http.scaladsl.model.ContentTypes
import better.files.File.OpenOptions
import cats.effect.IO
import com.google.api.gax.retrying.RetrySettings
import com.google.api.services.storage.StorageScopes
import com.google.api.services.storage.model.StorageObject
import com.google.auth.Credentials
import com.google.cloud.storage.Storage.BlobTargetOption
import com.google.cloud.storage.contrib.nio.{CloudStorageConfiguration, CloudStorageFileSystem, CloudStoragePath}
import com.google.cloud.storage.{Blob, BlobId, BlobInfo, StorageOptions}
import com.google.common.net.UrlEscapers
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.cloudsupport.gcp.gcs.GcsStorage
import cromwell.core.WorkflowOptions
import cromwell.core.path.{NioPath, Path, PathBuilder}
import cromwell.filesystems.gcs.GcsEnhancedRequest._
import cromwell.filesystems.gcs.GcsPathBuilder._
import cromwell.filesystems.gcs.GoogleUtil._
import mouse.all._
import org.apache.commons.codec.digest.DigestUtils

import java.io._
import java.net.URI
import scala.concurrent.{ExecutionContext, Future}
import scala.io.Codec
import scala.language.postfixOps
import scala.util.{Failure, Try}
object GcsPathBuilder {
  /*
   * Provides some level of validation of GCS bucket names
   * This is meant to alert the user early if they mistyped a gcs path in their workflow / inputs and not to validate
   * exact bucket syntax, which is done by GCS.
   * See https://cloud.google.com/storage/docs/naming for full spec
   */
  private val GcsBucketPattern =
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
    override def errorMessage: String =
      s"""
         |The path '$pathString' does not seem to be a valid GCS path.
         |Please check that it starts with gs:// and that the bucket and object follow GCS naming guidelines at
         |https://cloud.google.com/storage/docs/naming.
      """.stripMargin.replace("\n", " ").trim
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

  def validateGcsPath(string: String): GcsPathValidation =
    Try {
      val uri = URI.create(UrlEscapers.urlFragmentEscaper().escape(string))
      if (uri.getScheme == null) PossiblyValidRelativeGcsPath
      else if (uri.getScheme.equalsIgnoreCase(CloudStorageFileSystem.URI_SCHEME)) {
        if (uri.getHost == null) {
          softBucketParsing(string) map { ValidFullGcsPath(_, uri.getPath) } getOrElse InvalidFullGcsPath(string)
        } else ValidFullGcsPath(uri.getHost, uri.getPath)
      } else InvalidScheme(string)
    } recover { case t => UnparseableGcsPath(string, t) } get

  def isGcsPath(nioPath: NioPath): Boolean =
    nioPath.getFileSystem.provider().getScheme.equalsIgnoreCase(CloudStorageFileSystem.URI_SCHEME)

  def fromAuthMode(authMode: GoogleAuthMode,
                   applicationName: String,
                   retrySettings: RetrySettings,
                   cloudStorageConfiguration: CloudStorageConfiguration,
                   options: WorkflowOptions,
                   defaultProject: Option[String]
  )(implicit as: ActorSystem, ec: ExecutionContext): Future[GcsPathBuilder] =
    authMode.retryCredentials(options, List(StorageScopes.DEVSTORAGE_FULL_CONTROL)) map { credentials =>
      fromCredentials(credentials, applicationName, retrySettings, cloudStorageConfiguration, options, defaultProject)
    }

  def fromCredentials(credentials: Credentials,
                      applicationName: String,
                      retrySettings: RetrySettings,
                      cloudStorageConfiguration: CloudStorageConfiguration,
                      options: WorkflowOptions,
                      defaultProject: Option[String]
  ): GcsPathBuilder = {
    // Grab the google project from Workflow Options if specified and set
    // that to be the project used by the StorageOptions Builder. If it's not
    // specified use the default project mentioned in config file
    val project: Option[String] = options.get("google_project").toOption match {
      case Some(googleProject) => Option(googleProject)
      case None => defaultProject
    }

    val storageOptions = GcsStorage.gcsStorageOptions(credentials, retrySettings, project)

    // Create a com.google.api.services.storage.Storage
    // This is the underlying api used by com.google.cloud.storage
    // By bypassing com.google.cloud.storage, we can create low level requests that can be batched
    val apiStorage = GcsStorage.gcsStorage(applicationName, storageOptions)

    new GcsPathBuilder(apiStorage, cloudStorageConfiguration, storageOptions)
  }
}

class GcsPathBuilder(apiStorage: com.google.api.services.storage.Storage,
                     cloudStorageConfiguration: CloudStorageConfiguration,
                     storageOptions: StorageOptions
) extends PathBuilder {
  private[gcs] val projectId = storageOptions.getProjectId
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
  def build(string: String): Try[GcsPath] =
    validateGcsPath(string) match {
      case ValidFullGcsPath(bucket, path) =>
        Try {
          val fileSystem = CloudStorageFileSystem.forBucket(bucket, cloudStorageConfiguration, storageOptions)
          val cloudStoragePath = fileSystem.getPath(path)
          GcsPath(cloudStoragePath, apiStorage, cloudStorage, projectId)
        }
      case PossiblyValidRelativeGcsPath =>
        Failure(new IllegalArgumentException(s"""Path "$string" does not have a gcs scheme"""))
      case invalid: InvalidGcsPath => Failure(new IllegalArgumentException(invalid.errorMessage))
    }

  override def name: String = "Google Cloud Storage"
}

case class GcsPath private[gcs] (nioPath: NioPath,
                                 apiStorage: com.google.api.services.storage.Storage,
                                 cloudStorage: com.google.cloud.storage.Storage,
                                 projectId: String
) extends Path {

  val filesystemTypeKey = "gcs"

  lazy val objectBlobId: Try[BlobId] = Try {
    val bucketName = cloudStoragePath.bucket
    val objectName = cloudStoragePath.toRealPath().toString
    if (!Option(objectName).exists(_.trim.nonEmpty)) {
      throw new IllegalArgumentException(
        /*
        Instead of returning a 'Not Found' error, calling the Google Cloud Storage API without an object lists
        all objects in the bucket: https://cloud.google.com/storage/docs/listing-objects"

        List all objects:
        curl -X GET -H "Authorization: Bearer $(gcloud auth print-access-token $USER@broadinstitute.org)" \
          "https://storage.googleapis.com/storage/v1/b/hsd_salmon_index/o"

        List all objects, still:
        curl -X GET -H "Authorization: Bearer $(gcloud auth print-access-token $USER@broadinstitute.org)" \
          "https://storage.googleapis.com/storage/v1/b/hsd_salmon_index/o/"

        List a pseudo-directory named '/', and thus returns Not Found:
        curl -X GET -H "Authorization: Bearer $(gcloud auth print-access-token $USER@broadinstitute.org)" \
          "https://storage.googleapis.com/storage/v1/b/hsd_salmon_index/o//"

        See PROD-444 and linked tickets for more background.
         */
        s"GcsPath '$this' is bucket only and does not specify an object blob."
      )
    }
    BlobId.of(bucketName, objectName)
  }

  lazy val bucketOrObjectBlobId: Try[BlobId] =
    Try {
      val bucketName = cloudStoragePath.bucket
      val objectName = cloudStoragePath.toRealPath().toString
      BlobId.of(bucketName, objectName)
    }

  override protected def newPath(nioPath: NioPath): GcsPath = GcsPath(nioPath, apiStorage, cloudStorage, projectId)

  override def pathAsString: String = {
    val host = cloudStoragePath.bucket().stripSuffix("/")
    val path = cloudStoragePath.toString.stripPrefix("/")
    s"${CloudStorageFileSystem.URI_SCHEME}://$host/$path"
  }

  override def writeContent(content: String)(openOptions: OpenOptions, codec: Codec, compressPayload: Boolean)(implicit
    ec: ExecutionContext
  ): GcsPath.this.type = {
    def request(withProject: Boolean) = {
      val builder = BlobInfo.newBuilder(objectBlobId.get).setContentType(ContentTypes.`text/plain(UTF-8)`.value)
      if (compressPayload) {
        builder.setContentEncoding("gzip")
      }

      val contentByteArray = content.getBytes(codec.charSet)
      cloudStorage.create(
        builder.build(),
        if (compressPayload) gzipByteArray(contentByteArray) else contentByteArray,
        withProject.option(userProjectBlobTarget).toList.flatten: _*
      )
      this
    }

    // Since the NIO interface is synchronous we need to run this synchronously here. It is however wrapped in a Future
    // in the NioFlow so we don't need to worry about exceptions
    runOnEc(recoverFromProjectNotProvided(this, request)).unsafeRunSync()
    this
  }

  override def mediaInputStream(implicit ec: ExecutionContext): InputStream = {
    def request(withProject: Boolean) =
      // Use apiStorage here instead of cloudStorage, because apiStorage throws now if the bucket has requester pays,
      // whereas cloudStorage creates the input stream anyway and only throws one `read` is called (which only happens in NioFlow)
      apiStorage
        .objects()
        .get(objectBlobId.get.getBucket, objectBlobId.get.getName)
        .setUserProject(withProject.option(projectId).orNull)
        .executeMediaAsInputStream()
    // Since the NIO interface is synchronous we need to run this synchronously here. It is however wrapped in a Future
    // in the NioFlow so we don't need to worry about exceptions
    runOnEc(recoverFromProjectNotProvided(this, request)).unsafeRunSync()
  }

  override def pathWithoutScheme: String = {
    val gcsPath = cloudStoragePath
    gcsPath.bucket + gcsPath.toAbsolutePath.toString
  }

  def cloudStoragePath: CloudStoragePath = nioPath match {
    case gcsPath: CloudStoragePath => gcsPath
    case _ => throw new RuntimeException(s"Internal path was not a cloud storage path: $nioPath")
  }

  private def runOnEc[A](io: IO[A])(implicit ec: ExecutionContext): IO[A] = for {
    _ <- IO.shift(ec)
    result <- io
  } yield result

  private def userProjectBlobTarget: List[BlobTargetOption] = List(BlobTargetOption.userProject(projectId))
}

object GcsPath {

  // When generating "identity" hashes for GCS files (see cromwell.core.callcaching.FileHashStrategy), we use the
  // blob id, which includes bucket name, blob name, and generation. This slash-delimited string is md5-hashed to
  // obtain a string of a consistent length. This is used for call caching. Logic for generating this fingerprint
  // is centralized here because different codepaths are dealing with different representations of the GS object.

  def getBlobFingerprint(blob: StorageObject): Option[String] =
    makeFingerprint[StorageObject](
      blob,
      b => Option(b.getBucket),
      b => Option(b.getName),
      b =>
        // Prevent "Suspicious application of an implicit view" compiler warning by
        // explicitly converting java.lang.Long to scala Long
        Option(Long.unbox(b.getGeneration))
    )

  def getBlobFingerprint(blob: Blob): Option[String] =
    makeFingerprint[Blob](
      blob,
      b => Option(b.getBucket),
      b => Option(b.getName),
      b =>
        // Prevent "Suspicious application of an implicit view" compiler warning by
        // explicitly converting java.lang.Long to scala Long
        Option(Long.unbox(b.getGeneration))
    )

  private def makeFingerprint[T](blob: T,
                                 getBucket: T => Option[String],
                                 getName: T => Option[String],
                                 getGeneration: T => Option[Long]
  ): Option[String] =
    for {
      bucket <- getBucket(blob)
      name <- getName(blob)
      generation <- getGeneration(blob)
      id = List(bucket, name, generation.toString).mkString("/")
      _ = println(id)
    } yield DigestUtils.md5Hex(id).toUpperCase
}
