package cromwell.filesystems.s3

import java.net.URI

import com.google.common.net.UrlEscapers
import cromwell.core.path.PathBuilder
import cromwell.filesystems.s3.S3PathBuilder._
import org.lerch.s3fs.{AmazonS3ClientFactory, S3FileSystemProvider}
import scala.language.postfixOps
import scala.util.{Failure, Try}

object GenericS3PathBuilder {

  val S3BucketPattern =
    """
      (?x)                                      # Turn on comments and whitespace insensitivity
      ^s3://
      (                                         # Begin capturing group for bucket name
        [a-z0-9][a-z0-9-_\.]+[a-z0-9]           # Regex for bucket name - soft validation, see comment above
      )                                         # End capturing group for bucket name
      (?:
        /.*                                     # No validation here
      )?
    """.trim.r

  sealed trait S3PathValidation
  case class ValidFullS3Path(bucket: String, path: String) extends S3PathValidation
  case object PossiblyValidRelativeS3Path extends S3PathValidation
  sealed trait InvalidS3Path extends S3PathValidation {
    def pathString: String
    def errorMessage: String
  }
  final case class InvalidScheme(pathString: String) extends InvalidS3Path {
    override def errorMessage: String = s"S3 URIs must have 's3' scheme: $pathString"
  }
  final case class InvalidFullS3Path(pathString: String) extends InvalidS3Path {
    override def errorMessage: String = {
      s"""
         |The path '$pathString' does not seem to be a valid S3 path.
         |Please check that it starts with s3:// and that the bucket and object follow S3 naming guidelines at
         |https://docs.aws.amazon.com/AmazonS3/latest/dev/BucketRestrictions.html
      """.stripMargin.replace("\n", " ").trim
    }
  }
  final case class UnparseableS3Path(pathString: String, throwable: Throwable) extends InvalidS3Path {
    override def errorMessage: String =
      List(s"The specified S3 path '$pathString' does not parse as a URI.", throwable.getMessage).mkString("\n")
  }

  // Tries to extract a bucket name out of the provided string
  private def softBucketParsing(string: String): Option[String] = string match {
    case S3BucketPattern(bucket) => Option(bucket)
    case _ => None
  }

  def pathToUri(string: String): URI =
    URI.create(UrlEscapers.urlFragmentEscaper.escape(string))

  def validatePath(string: String): S3PathValidation = {
    Try {
      val uri = pathToUri(string)
      if (uri.getScheme == null) { PossiblyValidRelativeS3Path }
      else if (uri.getScheme.equalsIgnoreCase("s3")) {
        if (uri.getHost == null) {
          softBucketParsing(string) map { ValidFullS3Path(_, uri.getPath) } getOrElse InvalidFullS3Path(string)
        } else { ValidFullS3Path(uri.getHost, uri.getPath) }
      } else { InvalidScheme(string) }
    } recover { case t => UnparseableS3Path(string, t) } get
  }

  /*
  def fromAuthMode(authMode: AwsAuthMode,
                   configuration: S3Configuration,
                   options: WorkflowOptions,
                   storageRegion: Option[Region])(implicit ec: ExecutionContext): Future[GenericS3PathBuilder] = {
    val provider = authMode.provider()

    // Other backends needed retry here. In case we need retry, we'll return
    // a future. This will allow us to add capability without changing signature
    Future(fromProvider(provider,
      configuration,
      options,
      storageRegion
    ))
  }
  */

  def fromStaticMode(endpoint: String, accessKey: String, secretKey: String): GenericS3PathBuilder = {
    new GenericS3PathBuilder(endpoint, accessKey, secretKey)
  }

    /*
  def fromProvider(provider: AwsCredentialsProvider,
                   configuration: S3Configuration,
                   options: WorkflowOptions,
                   storageRegion: Option[Region]): GenericS3PathBuilder = {
    new GenericS3PathBuilder(configuration)
  }
  */
}

class GenericS3PathBuilder(endpoint: String, accessKey: String, secretKey: String) extends PathBuilder {
  // Tries to create a new S3Path from a String representing an absolute s3 path: s3://<bucket>[/<key>].
  def build(string: String): Try[S3Path] = {
    validatePath(string) match {
      case ValidFullS3Path(bucket, path) =>
        Try {
          val s3Path = new S3FileSystemProvider()
            .getFileSystem(URI.create("s3:////"), System.getenv)
            .getPath(s"""/$bucket/$path""")
          S3Path(s3Path, bucket,
            new AmazonS3ClientFactory().getS3Client(URI.create("s3:////"), System.getProperties))
        }
      case PossiblyValidRelativeS3Path => Failure(new IllegalArgumentException(s"$string does not have a s3 scheme"))
      case invalid: InvalidS3Path => Failure(new IllegalArgumentException(invalid.errorMessage))
    }
  }

  override def name: String = "s3"
}
