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
package cromwell.filesystems.aws

import java.net.URI
import java.nio.file.Paths
import akka.actor.ActorSystem
import com.google.common.net.UrlEscapers
import software.amazon.awssdk.core.auth.AwsCredentials
import software.amazon.awssdk.services.s3.{S3AdvancedConfiguration,S3Client}
import cromwell.cloudsupport.aws.auth.AwsAuthMode
import cromwell.cloudsupport.aws.s3.S3Storage
import cromwell.core.WorkflowOptions
import cromwell.core.path.{NioPath, Path, PathBuilder}
import cromwell.filesystems.aws.S3PathBuilder._

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

  def validatePath(string: String): S3PathValidation = {
    Try {
      val uri = URI.create(UrlEscapers.urlFragmentEscaper.escape(string))
      if (uri.getScheme == null) { PossiblyValidRelativeS3Path }
      else if (uri.getScheme.equalsIgnoreCase("s3")) {
        if (uri.getHost == null) {
          softBucketParsing(string) map { ValidFullS3Path(_, uri.getPath) } getOrElse InvalidFullS3Path(string)
        } else { ValidFullS3Path(uri.getHost, uri.getPath) }
      } else { InvalidScheme(string) }
    } recover { case t => UnparseableS3Path(string, t) } get
  }

  // TODO: Evaluate and remove or re-implement
  // We aren't using a S3 Filesystem provider. Leaving this commented because
  // someone might be dependent on it later
  // def isS3Path(nioPath: NioPath): Boolean = {
  //   nioPath.getFileSystem.provider().getScheme.equalsIgnoreCase("s3")
  // }

  def fromAuthMode(authMode: AwsAuthMode,
                   configuration: S3AdvancedConfiguration,
                   options: WorkflowOptions)(implicit ec: ExecutionContext): Future[S3PathBuilder] = {
    val credentials = authMode.credential((key: String) => options.get(key).get)

    // Other backends needed retry here. In case we need retry, we'll return
    // a future. This will allow us to add capability without changing signature
    Future(fromCredentials(credentials,
      configuration,
      options
    ))
  }

  def fromCredentials(credentials: AwsCredentials,
                      configuration: S3AdvancedConfiguration,
                      options: WorkflowOptions): S3PathBuilder = {
    new S3PathBuilder(S3Storage.s3Client(credentials), configuration)
  }
}

class S3PathBuilder(client: S3Client,
                     configuration: S3AdvancedConfiguration
                     ) extends PathBuilder {


  // Tries to create a new S3Path from a String representing an absolute s3 path: s3://<bucket>[/<key>].
  def build(string: String): Try[S3Path] = {
    validatePath(string) match {
      case ValidFullS3Path(bucket, path) =>
        Try {
          // TODO: Verify that the s3:// doesn't get in the way of this...
          S3Path(Paths.get(path), bucket, client)
        }
      case PossiblyValidRelativeS3Path => Failure(new IllegalArgumentException(s"$string does not have a s3 scheme"))
      case invalid: InvalidS3Path => Failure(new IllegalArgumentException(invalid.errorMessage))
    }
  }

  override def name: String = "S3"
}

case class S3Path private[aws](nioPath: NioPath,
                               bucket: String,
                               client: S3Client,
                               ) extends Path {
  override protected def newPath(nioPath: NioPath): S3Path = S3Path(nioPath, bucket, client)

  override def pathAsString: String = {
    val host = bucket.stripSuffix("/")
    val path = nioPath.toString.stripPrefix("/")
    s"s3://$host/$path"
  }

  override def pathWithoutScheme: String = {
    bucket + nioPath.toAbsolutePath.toString
  }

  def key: String = {
    nioPath.toAbsolutePath.toString
  }
}
