package cromwell.filesystems.s3

import java.net.URI
import java.util.Properties

import com.google.common.net.UrlEscapers

import cromwell.filesystems.s3.S3PathBuilder._
import org.lerch.s3fs.{AmazonS3Factory,S3FileSystemProvider}
import scala.language.postfixOps
import scala.util.{Failure, Try}

import cromwell.core.path.{NioPath, Path, PathBuilder}
import org.lerch.s3fs.util.S3Utils

import software.amazon.awssdk.services.s3.S3Client


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

  def fromStaticMode(endpoint: String, accessKey: String, secretKey: String): GenericS3PathBuilder = {
    new GenericS3PathBuilder(endpoint, accessKey, secretKey)
  }

}

class GenericS3PathBuilder(endpoint: String, accessKey: String, secretKey: String) extends PathBuilder {
  // Tries to create a new S3Path from a String representing an absolute s3 path: s3://<bucket>[/<key>].
  def build(string: String): Try[GenericS3Path] = {
    validatePath(string) match {
      case ValidFullS3Path(bucket, path) =>
        Try {
          val props = new Properties()
          props.setProperty( AmazonS3Factory.ACCESS_KEY, accessKey)
          props.setProperty( AmazonS3Factory.SECRET_KEY, secretKey)
          props.setProperty( AmazonS3Factory.REGION, "us-west")
          val s3Path = new S3FileSystemProvider()
            .createFileSystem(URI.create(endpoint), props)
            .getPath(s"""/$bucket/$path""")
          val client = s3Path.getFileSystem().getClient()
          GenericS3Path(s3Path, bucket, client)
        }
      case PossiblyValidRelativeS3Path => Failure(new IllegalArgumentException(s"$string does not have a s3 scheme"))
      case invalid: InvalidS3Path => Failure(new IllegalArgumentException(invalid.errorMessage))
    }
  }

  override def name: String = "s3"
}



case class GenericS3Path private[s3](nioPath: NioPath,
                               bucket: String,
                               client: S3Client
                               ) extends Path {
  override protected def newPath(nioPath: NioPath): GenericS3Path = GenericS3Path(nioPath, bucket, client)

  override def pathAsString: String = s"s3://$pathWithoutScheme"

  override def pathWithoutScheme: String = {
    (new URI(safeAbsolutePath)).getPath().stripPrefix("/")
  }
  def key: String = safeAbsolutePath

  /*
    This is a bit of a headache. We could make S3Path in S3PathBuilder as that can fulfill all usages of
    nioPath. However cromwell.core.path.Path requires the NioPath argument. Further we already know that
    nioPath is a valid S3Path by nature of how it was created in S3PathBuilder.build. However, better safe than
    sorry
   */
  lazy val s3Path = nioPath match {
    case s3Path: org.lerch.s3fs.S3Path => s3Path
    case _ => throw new RuntimeException("Internal path was not an S3 path: " + nioPath)
  }

  lazy val eTag = new S3Utils().getS3ObjectSummary(s3Path).eTag()

  /** Gets an absolute path for multiple forms of input. The FS provider does
   *  not support "toAbsolutePath" on forms such as "mypath/" or "foo.bar"
   *  So this function will prepend a forward slash for input that looks like this
   *  while leaving properly rooted input or input beginning with s3:// alone
   */
  def safeAbsolutePath: String = {
    val originalPath = s3Path.toString
    if (originalPath.startsWith("s3")) return s3Path.toAbsolutePath.toString
    originalPath.charAt(0) match {
      case '/' =>  s3Path.toAbsolutePath.toString
      case _ => s3Path.resolve(s"/$bucket/$originalPath").toAbsolutePath.toString
    }
  }
}
