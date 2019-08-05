/*
 * Copyright 2018 Amazon.com, Inc. or its affiliates.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 *  BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 *  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 *  THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 *  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 *  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 *  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */
package cromwell.filesystems.s3

import java.net.URI

import com.google.common.net.UrlEscapers
import software.amazon.awssdk.auth.credentials.AwsCredentials
import software.amazon.awssdk.services.s3.{S3Client, S3Configuration}
import cromwell.cloudsupport.aws.auth.AwsAuthMode
import cromwell.cloudsupport.aws.s3.S3Storage
import cromwell.core.WorkflowOptions
import cromwell.core.path.{NioPath, Path, PathBuilder}
import cromwell.filesystems.s3.S3PathBuilder._
import org.lerch.s3fs.S3FileSystemProvider
import org.lerch.s3fs.util.S3Utils
import software.amazon.awssdk.regions.Region

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Try}

object S3PathBuilder {

  // Provides some level of validation of bucket names
  // This is meant to alert the user early if they mistyped a path in their workflow / inputs and not to validate
  // exact bucket syntax.
  // See https://docs.aws.amazon.com/AmazonS3/latest/dev/BucketRestrictions.html
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

  def fromAuthMode(authMode: AwsAuthMode,
                   configuration: S3Configuration,
                   options: WorkflowOptions,
                   storageRegion: Option[Region])(implicit ec: ExecutionContext): Future[S3PathBuilder] = {
    val credentials = authMode.credential((key: String) => options.get(key).get)

    // Other backends needed retry here. In case we need retry, we'll return
    // a future. This will allow us to add capability without changing signature
    Future(fromCredentials(credentials,
      configuration,
      options,
      storageRegion
    ))
  }

  def fromCredentials(credentials: AwsCredentials,
                      configuration: S3Configuration,
                      options: WorkflowOptions,
                      storageRegion: Option[Region]): S3PathBuilder = {
    new S3PathBuilder(S3Storage.s3Client(credentials, storageRegion), configuration)
  }
}

class S3PathBuilder(client: S3Client,
                     configuration: S3Configuration
                     ) extends PathBuilder {
  // Tries to create a new S3Path from a String representing an absolute s3 path: s3://<bucket>[/<key>].
  def build(string: String): Try[S3Path] = {
    validatePath(string) match {
      case ValidFullS3Path(bucket, path) =>
        Try {
          // TODO: System.getenv needs to turn into a full Auth thingy
          // TODO: This assumes the "global endpoint". Need to handle other endpoints
          val s3Path = new S3FileSystemProvider()
            .getFileSystem(URI.create("s3:////"), System.getenv, client)
            .getPath(s"""/$bucket/$path""")
          S3Path(s3Path, bucket, client)
        }
      case PossiblyValidRelativeS3Path => Failure(new IllegalArgumentException(s"$string does not have a s3 scheme"))
      case invalid: InvalidS3Path => Failure(new IllegalArgumentException(invalid.errorMessage))
    }
  }

  override def name: String = "s3"
}

case class S3Path private[s3](nioPath: NioPath,
                               bucket: String,
                               client: S3Client
                               ) extends Path {
  override protected def newPath(nioPath: NioPath): S3Path = S3Path(nioPath, bucket, client)

  override def pathAsString: String = s"s3://$pathWithoutScheme"

  override def pathWithoutScheme: String = safeAbsolutePath.stripPrefix("s3://s3.amazonaws.com/")

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
