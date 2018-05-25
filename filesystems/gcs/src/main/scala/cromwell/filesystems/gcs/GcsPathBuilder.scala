package cromwell.filesystems.gcs

import java.net.URI

import akka.actor.ActorSystem
import com.google.api.gax.retrying.RetrySettings
import com.google.auth.Credentials
import com.google.cloud.storage.contrib.nio.{CloudStorageConfiguration, CloudStorageFileSystem, CloudStoragePath}
import com.google.cloud.storage.{BlobId, StorageOptions}
import com.google.common.net.UrlEscapers
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.cloudsupport.gcp.gcs.GcsStorage
import cromwell.core.WorkflowOptions
import cromwell.core.path.{NioPath, Path, PathBuilder}
import cromwell.filesystems.gcs.GcsPathBuilder._
import cromwell.filesystems.gcs.GoogleUtil._

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Try}

object GcsPathBuilder {
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
                   defaultProject: Option[String])(implicit as: ActorSystem, ec: ExecutionContext): Future[GcsPathBuilder] = {
    authMode.retryCredential(options) map { credentials =>
      fromCredentials(credentials,
        applicationName,
        retrySettings,
        cloudStorageConfiguration,
        options,
        defaultProject
      )
    }
  }

  def fromCredentials(credentials: Credentials,
                      applicationName: String,
                      retrySettings: RetrySettings,
                      cloudStorageConfiguration: CloudStorageConfiguration,
                      options: WorkflowOptions,
                      defaultProject: Option[String]): GcsPathBuilder = {
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

    new GcsPathBuilder(apiStorage, cloudStorageConfiguration, storageOptions)
  }
}

class GcsPathBuilder(apiStorage: com.google.api.services.storage.Storage,
                     cloudStorageConfiguration: CloudStorageConfiguration,
                     storageOptions: StorageOptions) extends PathBuilder {

  private[gcs] val projectId = storageOptions.getProjectId

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
        Try {
          val fileSystem = CloudStorageFileSystem.forBucket(bucket, cloudStorageConfiguration, storageOptions)
          val cloudStoragePath = fileSystem.getPath(path)
          val cloudStorage = storageOptions.getService
          GcsPath(cloudStoragePath, apiStorage, cloudStorage, projectId)
        }
      case PossiblyValidRelativeGcsPath => Failure(new IllegalArgumentException(s"$string does not have a gcs scheme"))
      case invalid: InvalidGcsPath => Failure(new IllegalArgumentException(invalid.errorMessage))
    }
  }

  override def name: String = "Google Cloud Storage"
}

case class GcsPath private[gcs](nioPath: NioPath,
                                apiStorage: com.google.api.services.storage.Storage,
                                cloudStorage: com.google.cloud.storage.Storage,
                                projectId: String) extends Path {
  lazy val blob = BlobId.of(cloudStoragePath.bucket, cloudStoragePath.toRealPath().toString)

  override protected def newPath(nioPath: NioPath): GcsPath = GcsPath(nioPath, apiStorage, cloudStorage, projectId)

  override def pathAsString: String = {
    val host = cloudStoragePath.bucket().stripSuffix("/")
    val path = cloudStoragePath.toString.stripPrefix("/")
    s"${CloudStorageFileSystem.URI_SCHEME}://$host/$path"
  }

  /*
   * This is commented out for now as it forces service accounts and users to have a specific role with permssion
   * to use Cromwell. However it is needed to enable Requester pays. This code can be uncommented / adjusted when fully
   * enable requester pays.
  
  /***
    * This method needs to be overridden to make it work with requester pays. We need to go around Nio
    * as currently it doesn't support to set the billing project id. Google Cloud Storage already has
    * code in place to set the billing project (inside HttpStorageRpc) but Nio does not pass it even though
    * it's available. In future when it is supported, remove this method and wire billing project into code
    * wherever necessary
    */
  override def mediaInputStream: InputStream = {
    Try{
      apiStorage.objects().get(blob.getBucket, blob.getName).setUserProject(projectId).executeMediaAsInputStream()
    } match {
      case Success(inputStream) => inputStream
      case Failure(e) => throw new IOException(s"Failed to open an input stream for $pathAsString", e)
    }
  }

  /***
    * This method needs to be overridden to make it work with requester pays. We need to go around Nio
    * as currently it doesn't support to set the billing project id. Google Cloud Storage already has
    * code in place to set the billing project (inside HttpStorageRpc) but Nio does not pass it even though
    * it's available. In future when it is supported, remove this method and wire billing project into code
    * wherever necessary
    */
  override def readContentAsString(implicit codec: Codec): String = {
    openInputStream(readLinesAsString)
  }

  private def readLinesAsString(inputStream: InputStream)(implicit codec: Codec): String = {
    val byteArray = Stream.continually(inputStream.read).takeWhile(_ != -1).map(_.toByte).toArray
    new String(byteArray, Charset.forName(codec.name))
  }

  /***
    * This method needs to be overridden to make it work with requester pays. We need to go around Nio
    * as currently it doesn't support to set the billing project id. Google Cloud Storage already has
    * code in place to set the billing project (inside HttpStorageRpc) but Nio does not pass it even though
    * it's available. In future when it is supported, remove this method and wire billing project into code
    * wherever necessary
    */
  override def readAllLinesInFile(implicit codec: Codec): Traversable[String] = {
    openInputStream(readLinesAsStringTraversable)
  }

  private def readLinesAsStringTraversable(inputStream: InputStream)(implicit codec: Codec): Traversable[String] = {
    val reader = new BufferedReader(new InputStreamReader(inputStream, codec.name))
    Stream.continually(reader.readLine()).takeWhile(_ != null).toList
  }

  private def openInputStream[A](f: InputStream => A): A = {
    val storageObject = apiStorage.objects().get(blob.getBucket, blob.getName).setUserProject(projectId)

    val output = tryWithResource(() => storageObject.executeMediaAsInputStream())(inputStream => f(inputStream))

    output.getOrElse(throw new IOException(s"Failed to open an input stream for $pathAsString"))
  }
  */

  override def pathWithoutScheme: String = {
    val gcsPath = cloudStoragePath
    gcsPath.bucket + gcsPath.toAbsolutePath.toString
  }

  def cloudStoragePath: CloudStoragePath = nioPath match {
    case gcsPath: CloudStoragePath => gcsPath
    case _ => throw new RuntimeException(s"Internal path was not a cloud storage path: $nioPath")
  }
}
