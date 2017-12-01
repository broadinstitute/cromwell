package cromwell.filesystems.oss

import java.io.{BufferedReader, ByteArrayInputStream, FileInputStream, InputStreamReader}
import java.net.URI

import com.aliyun.oss.OSSClient
import com.google.cloud.storage.contrib.nio.{CloudStorageConfiguration, CloudStorageFileSystem, CloudStoragePath}
import com.google.common.net.UrlEscapers
import cromwell.core.WorkflowOptions
import cromwell.core.path.BetterFileMethods.{LinkOptions, OpenOptions}
import cromwell.core.path.{NioPath, Path, PathBuilder}
import cromwell.filesystems.oss.OssPathBuilder._

import scala.io.Codec
import scala.language.postfixOps
import scala.util.control.Breaks
import scala.util.matching.Regex
import scala.util.{Failure, Try}

object OssPathBuilder {
  /* 
    * Provides some level of validation of OSS bucket names
    * This is meant to alert the user early if they mistyped a gcs path in their workflow / inputs and not to validate
    * exact bucket syntax, which is done by OSS.
  */
  val DefaultCloudStorageConfiguration = {
    CloudStorageConfiguration.builder()
      .blockSize(1024)
      .permitEmptyPathComponents(true)
      .stripPrefixSlash(true)
      .usePseudoDirectories(true)
      .build()
  }

  val URI_SCHEME = "oss"
  val OssBucketPattern:Regex =
    """
      (?x)                                      # Turn on comments and whitespace insensitivity
      ^oss://
      (                                         # Begin capturing group for gcs bucket name
        [a-z0-9][a-z0-9-_\\.]+[a-z0-9]          # Regex for bucket name - soft validation, see comment above
      )                                         # End capturing group for gcs bucket name
      (?:
        /.*                                     # No validation here
      )?
    """.trim.r

  sealed trait OssPathValidation
  case class ValidFullOssPath(bucket: String, path: String) extends OssPathValidation
  case object PossiblyValidRelativeOssPath extends OssPathValidation
  sealed trait InvalidOssPath extends OssPathValidation {
    def pathString: String
    def errorMessage: String
  }

  final case class InvalidScheme(pathString: String) extends InvalidOssPath {
    def errorMessage = s"OSS URIs must have 'oss' scheme: $pathString"
  }
  final case class InvalidFullGcsPath(pathString: String) extends InvalidOssPath {
    def errorMessage = {
      s"""
         |The path '$pathString' does not seem to be a valid OSS path.
         |Please check that it starts with oss:// and that the bucket and object follow OSS naming guidelines.
      """.stripMargin.replaceAll("\n", " ").trim
    }
  }
  final case class UnparseableOssPath(pathString: String, throwable: Throwable) extends InvalidOssPath {
    def errorMessage: String =
      List(s"The specified OSS path '$pathString' does not parse as a URI.", throwable.getMessage).mkString("\n")
  }

  private def softBucketParsing(string: String): Option[String] = string match {
    case OssBucketPattern(bucket) => Option(bucket)
    case _ => None
  }

  def validateOssPath(string: String): OssPathValidation = {
    Try {
      val uri = URI.create(UrlEscapers.urlFragmentEscaper().escape(string))
      if (uri.getScheme == null) PossiblyValidRelativeOssPath
      else if (uri.getScheme.equalsIgnoreCase(URI_SCHEME)) {
        if (uri.getHost == null) {
          softBucketParsing(string) map { ValidFullOssPath(_, uri.getPath) } getOrElse InvalidFullGcsPath(string)
        } else ValidFullOssPath(uri.getHost, uri.getPath)
      } else InvalidScheme(string)
    } recover { case t => UnparseableOssPath(string, t) } get
  }

  def isOssPath(nioPath: NioPath): Boolean = {
    nioPath.getFileSystem.provider().getScheme.equalsIgnoreCase(URI_SCHEME)
  }

  def fromCredentials(endpoint: String,
                      accessId: String,
                      accessKey: String,
                      securityToken: Option[String],
                      options: WorkflowOptions): OssPathBuilder = {
    val ossClient: OSSClient = securityToken match {
      case Some(token) => new OSSClient(endpoint, accessId, accessKey, token)
      case None => new OSSClient(endpoint, accessId, accessKey)
    }

    new OssPathBuilder(ossClient)
  }
}

class OssPathBuilder(ossClient: OSSClient) extends PathBuilder {

  val client = ossClient
  def build(string: String): Try[OssPath] = {
    validateOssPath(string) match {
      case ValidFullOssPath(bucket, path) =>
        Try {
          val fileSystem = CloudStorageFileSystem.forBucket(bucket, DefaultCloudStorageConfiguration)
          val cloudStoragePath = fileSystem.getPath(path)
          new OssPath(cloudStoragePath, client)
        }
      case PossiblyValidRelativeOssPath => Failure(new IllegalArgumentException(s"$string does not have a oss scheme"))
      case invalid: InvalidOssPath => Failure(new IllegalArgumentException(invalid.errorMessage))
    }
  }

  override def name: String = "Oss"
}

case class BucketAndObj(bucket: String, obj: String)

case class OssPath private[oss](nioPath: NioPath,
                                ossClient: OSSClient,
                               ) extends Path {

  override protected def newPath(path: NioPath): OssPath = {
    OssPath(path, ossClient)
  }

  override def pathAsString: String = {
    val host = cloudStoragePath.bucket().stripSuffix("/")
    val path = cloudStoragePath.toString.stripPrefix("/")
    s"${OssPathBuilder.URI_SCHEME}://$host/$path"
  }

  override def pathWithoutScheme: String = {
    val gcsPath = cloudStoragePath
    gcsPath.bucket + gcsPath.toAbsolutePath.toString
  }

  def bucket: String = {
    return cloudStoragePath.bucket().stripSuffix("/")
  }

  def key: String = {
    return cloudStoragePath.toString.stripPrefix("/")
  }

  def cloudStoragePath: CloudStoragePath = nioPath match {
    case gcsPath: CloudStoragePath => gcsPath
    case _ => throw new RuntimeException(s"Internal path was not a cloud storage path: $nioPath")
  }

  override def size: Long = {
    val meta = ossClient.getObjectMetadata(bucket, key)
    meta.getContentLength
  }

  override def contentAsString(implicit codec: Codec): String = {
    val ossObject = ossClient.getObject(bucket, key)

    val reader = new BufferedReader(new InputStreamReader(ossObject.getObjectContent()))

    val ret = new StringBuilder
    val loop = new Breaks
    loop.breakable {
      while (true) {
        val line = reader.readLine()
        if (line == null) {
          loop.break
        }
        ret ++= (line + '\n')
      }
    }
    reader.close()
    ret.toString()
  }

  override def write(text: String)(implicit openOptions: OpenOptions = OpenOptions.default, codec: Codec): this.type = {
    openOptions match {
      case Seq(UploadFileOption) =>
        val inputStream = new FileInputStream(text)
        ossClient.putObject(bucket, key, inputStream)
        this
      case _ =>
        ossClient.putObject(bucket, key, new ByteArrayInputStream(text.getBytes()))
        this
    }
  }

  override def exists(implicit linkOptions: LinkOptions = LinkOptions.default): Boolean = {
    if (key.isEmpty)
      ossClient.doesBucketExist(bucket)
    else
      ossClient.doesObjectExist(bucket, key)
  }

  override def notExists(implicit linkOptions: LinkOptions = LinkOptions.default): Boolean = {
    !exists(linkOptions)
  }
}
