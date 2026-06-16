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
import cromwell.cloudsupport.aws.auth.AwsAuthMode
import cromwell.core.WorkflowOptions
import cromwell.core.path.{NioPath, Path, PathBuilder}
import cromwell.filesystems.s3.S3PathBuilder._
import cromwell.cloudsupport.aws.s3.S3Storage
import org.lerch.s3fs.{S3FileSystem, S3FileSystemProvider}
import org.lerch.s3fs.util.S3Utils
import software.amazon.awssdk.auth.credentials.AwsCredentialsProvider
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.s3.{S3Client, S3Configuration}

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
    override def errorMessage: String =
      s"""
         |The path '$pathString' does not seem to be a valid S3 path.
         |Please check that it starts with s3:// and that the bucket and object follow S3 naming guidelines at
         |https://docs.aws.amazon.com/AmazonS3/latest/dev/BucketRestrictions.html
      """.stripMargin.replace("\n", " ").trim
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

  def validatePath(string: String): S3PathValidation =
    Try {
      val uri = pathToUri(string)
      if (uri.getScheme == null) { PossiblyValidRelativeS3Path }
      else if (uri.getScheme.equalsIgnoreCase("s3")) {
        if (uri.getHost == null) {
          softBucketParsing(string) map { ValidFullS3Path(_, uri.getPath) } getOrElse InvalidFullS3Path(string)
        } else { ValidFullS3Path(uri.getHost, uri.getPath) }
      } else { InvalidScheme(string) }
    } recover { case t => UnparseableS3Path(string, t) } get

  /**
   * @param endpointUri optional custom S3-compatible endpoint URI (non-AWS, e.g. OVH, MinIO).
   */
  def fromAuthMode(authMode: AwsAuthMode,
                   configuration: S3Configuration,
                   options: WorkflowOptions,
                   storageRegion: Option[Region],
                   endpointUri: Option[URI] = None
  )(implicit ec: ExecutionContext): Future[S3PathBuilder] = {
    val provider = authMode.provider()
    // Other backends needed retry here. In case we need retry, we'll return
    // a future. This will allow us to add capability without changing signature
    Future(fromProvider(provider, configuration, options, storageRegion, endpointUri))
  }

  def fromProvider(provider: AwsCredentialsProvider,
                   configuration: S3Configuration,
                   options: WorkflowOptions,
                   storageRegion: Option[Region],
                   endpointUri: Option[URI] = None
  ): S3PathBuilder =
    new S3PathBuilder(provider, configuration, storageRegion, endpointUri)
}

/**
 * Builds S3Path instances for a given credentials provider, region and (optionally) a
 * custom S3-compatible endpoint.
 *
 * Previously this class held only `S3Configuration` and `build()` always fell back
 * to `System.getenv` / `System.getProperties` for credentials and endpoint, making it
 * impossible to use configured credentials or a non-AWS endpoint. The class now stores the
 * full auth context and injects it into both the s3fs-nio filesystem (for NIO path resolution)
 * and the AWS SDK v2 S3Client (for direct SDK operations stored on S3Path).
 *
 * @param provider     AWS SDK v2 credentials provider built from the configured auth mode
 * @param configuration S3Configuration (accelerate, dual-stack, path-style flags)
 * @param storageRegion optional region (AWS region or S3-compatible region label)
 * @param endpointUri  optional custom S3-compatible endpoint; when set, path-style access
 *                     is forced and the s3fs-nio filesystem is pointed at this host
 */
class S3PathBuilder(
    provider: AwsCredentialsProvider,
    configuration: S3Configuration,
    storageRegion: Option[Region],
    endpointUri: Option[URI]
) extends PathBuilder {

  /**
   * Lazily-initialised S3FileSystem backed by our own AWS SDK v2 S3Client.
   * Shared across all build() calls on this builder instance (one builder per workflow).
   *
   * s3fs-nio's AmazonS3Factory.getS3Client(URI, Properties) extracts the host from
   * the filesystem URI and then calls AWS SDK v2 builder.endpointOverride(uri) passing the
   * full URI — including the s3:// scheme — which the SDK rejects:
   *   "Custom endpoint 's3://...' was not a valid URI"
   * The SDK requires http:// or https:// for endpointOverride.
   *
   * S3FileSystemProvider also exposes a 3-arg overload
   *   createFileSystem(URI, Properties, S3Client)
   * that constructs S3FileSystem directly from a caller-supplied S3Client, completely
   * bypassing AmazonS3Factory. We pre-build the S3Client via S3Storage.s3Client() which
   * already calls endpointOverride() with the correct https:// URI.
   */
  private lazy val cachedFilesystem: S3FileSystem = {
    // Build our correctly-configured S3Client (endpointOverride with https:// URI,
    // pathStyleAccessEnabled forced when endpointUri is set).
    val s3Client = S3Storage.s3Client(configuration, provider, storageRegion, endpointUri)

    // Build a Properties map with our configured credentials for s3fs-nio's bookkeeping.
    // s3fs-nio uses s3fs_access_key / s3fs_secret_key as part of the filesystem cache key.
    // We do not need to copy System.getenv() here because the actual S3 client is already
    // pre-built (passed as the third argument to createFileSystem); AmazonS3Factory is never
    // invoked, so only the cache-key fields matter.
    val creds = provider.resolveCredentials()
    val props  = new java.util.Properties()
    props.put("s3fs_access_key", creds.accessKeyId())
    props.put("s3fs_secret_key", creds.secretAccessKey())

    // The URI is only used to derive the filesystem key; the actual endpoint is
    // driven by s3Client, not by this URI.
    val fsUri = endpointUri
      .map(ep => URI.create(s"s3://${ep.getHost}/"))
      .getOrElse(URI.create("s3:////"))

    new S3FileSystemProvider().createFileSystem(fsUri, props, s3Client)
  }

  // Tries to create a new S3Path from a String representing an absolute s3 path: s3://<bucket>[/<key>].
  def build(string: String): Try[S3Path] =
    validatePath(string) match {
      case ValidFullS3Path(bucket, path) =>
        Try {
          val nioPath = cachedFilesystem.getPath(s"""/$bucket/$path""")
          S3Path(nioPath, bucket, cachedFilesystem.getClient())
        }
      case PossiblyValidRelativeS3Path => Failure(new IllegalArgumentException(s"$string does not have a s3 scheme"))
      case invalid: InvalidS3Path => Failure(new IllegalArgumentException(invalid.errorMessage))
    }

  override def name: String = "s3"
}

case class S3Path private[s3] (nioPath: NioPath, bucket: String, client: S3Client) extends Path {
  val filesystemTypeKey = "s3"

  override protected def newPath(nioPath: NioPath): S3Path = S3Path(nioPath, bucket, client)

  override def pathAsString: String = s"s3://$pathWithoutScheme"

  // Previously, this was `stripPrefix("s3://s3.amazonaws.com/")` which broke for custom S3-compatible
  // endpoints (OVH, MinIO, etc.) where s3fs-nio encodes a different host in the path string.
  // Now strips s3://[any-host]/ generically so pathAsString always returns s3://bucket/key.
  override def pathWithoutScheme: String = safeAbsolutePath.replaceFirst("^s3://[^/]*/", "")

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
      case '/' => s3Path.toAbsolutePath.toString
      case _ => s3Path.resolve(s"/$bucket/$originalPath").toAbsolutePath.toString
    }
  }

  /**
   * Override createDirectories() as a no-op for S3 paths.
   *
   * S3 is a flat key-value namespace — "directories" are a purely virtual concept defined
   * by key prefixes. There is no API call required (or meaningful) to "create" a directory.
   * The base BetterFileMethods.createDirectories() delegates to better-files which calls the
   * s3fs-nio NIO FileSystemProvider.createDirectory(), which attempts to PUT a zero-byte
   * object with a trailing slash as a "directory marker". S3-compatible endpoints such as
   * OVH Ceph RadosGW reject these putObject calls with 403 Access Denied.
   *
   * Since directories don't physically exist in S3, creating them is unnecessary.
   * Cromwell's actual workflow files are written directly at their full key paths later;
   * those writes succeed without any parent "directory" having been pre-created.
   *
   * The subsequent addPermission() call in createPermissionedDirectories() will throw
   * IOException (S3 has no POSIX permission model), which is already caught and ignored
   * by EvenBetterPathMethods.createPermissionedDirectories().
   */
  override def createDirectories()(implicit attributes: better.files.File.Attributes = better.files.File.Attributes.default): this.type = {
    // S3 directories are virtual — no-op.
    this
  }
}
