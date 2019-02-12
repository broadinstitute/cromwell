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

import software.amazon.awssdk.auth.credentials.AnonymousCredentialsProvider
import cromwell.cloudsupport.aws.s3.S3Storage
import cromwell.core.Tags.AwsTest
import cromwell.core.path._
import cromwell.core.{TestKitSuite, WorkflowOptions}
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpecLike, Matchers}
import software.amazon.awssdk.regions.Region

class S3PathBuilderSpec extends TestKitSuite with FlatSpecLike with Matchers with PathBuilderSpecUtils {

  behavior of "S3PathBuilder"

  it should "build from valid credentials" taggedAs AwsTest in {
    // We aren't using workflow options...yet
    val wfOptionsWithProject = WorkflowOptions.empty

    S3PathBuilder.fromCredentials(
      AnonymousCredentialsProvider.create.resolveCredentials(),
      S3Storage.s3Configuration(),
      wfOptionsWithProject,
      Option(Region.US_EAST_1)
    )
  }

  it should behave like truncateCommonRoots(pathBuilder, pathsToTruncate, AwsTest)

  goodPaths foreach { goodPath =>
    it should behave like buildGoodPath(pathBuilder, goodPath, AwsTest)
  }

  badPaths foreach { badPath =>
    it should behave like buildBadPath(pathBuilder, badPath, AwsTest)
  }

  private def pathsToTruncate = Table(
    ("context", "file", "relative"),
    ("s3://bucket", "s3://bucket/path/to/file", "path/to/file"),
    ("s3://bucket/path/to/my/dir", "s3://bucket/path/to/my/dir/file", "file"),
    ("s3://bucket/path/to/my/dir", "s3://bucket/path/to/my/dir//file", "file"),
    // TODO: Investigate these next two. DefaultPathBuilder does something different from
    //       GcsPathBuilder. Need to determine what's the right test
    // ("s3://bucket/path/to/my//dir", "s3://bucket/path/to/my/dir/file", "dir/file"),
    // ("s3://bucket/path/to/my//dir", "s3://bucket/path/to/my/dir//file", "dir//file"),
    ("s3://bucket/path/to/my/dir", "s3://bucket/path/./to/my/dir/file", "./to/my/dir/file"),
    ("s3://bucket/path/to/my/dir/with/file", "s3://bucket/path/to/other/dir/with/file", "other/dir/with/file")
  )

  private def bucket = "mymadeupbucket"

  private def goodPaths = Seq(
    // GoodPath(
    //   description = "a path with spaces",
    //   path = s"s3://$bucket/hello/world/with spaces",
    //   normalize = false,
    //   pathAsString = s"s3://$bucket/hello/world/with spaces",
    //   pathWithoutScheme = s"$bucket/hello/world/with spaces",
    //   parent = s"s3://$bucket/hello/world/",
    //   getParent = s"s3://$bucket/hello/world/",
    //   root = s"s3://$bucket/",
    //   name = "with spaces",
    //   getFileName = s"s3://$bucket/with spaces",
    //   getNameCount = 3,
    //   isAbsolute = true),
    GoodPath(
      description = "a path with spaces",
      path = s"s3://$bucket/hello/world/with spaces",
      normalize = false,
      pathAsString = s"s3://$bucket/hello/world/with spaces",
      pathWithoutScheme = s"$bucket/hello/world/with spaces",
      parent = s"s3://$bucket/hello/world",
      getParent = s"s3://$bucket/hello/world",
      root = s"s3://$bucket/",
      name = "with spaces",
      getFileName = s"s3://$bucket/with spaces",
      getNameCount = 3,
      isAbsolute = true),

    // GoodPath(
    //   description = "a path with non-ascii",
    //   path = s"s3://$bucket/hello/world/with non ascii £€",
    //   normalize = false,
    //   pathAsString = s"s3://$bucket/hello/world/with non ascii £€",
    //   pathWithoutScheme = s"$bucket/hello/world/with non ascii £€",
    //   parent = s"s3://$bucket/hello/world/",
    //   getParent = s"s3://$bucket/hello/world/",
    //   root = s"s3://$bucket/",
    //   name = "with non ascii £€",
    //   getFileName = s"s3://$bucket/with non ascii £€",
    //   getNameCount = 3,
    //   isAbsolute = true),
    GoodPath(
      description = "a path with non-ascii",
      path = s"s3://$bucket/hello/world/with non ascii £€",
      normalize = false,
      pathAsString = s"s3://$bucket/hello/world/with non ascii £€",
      pathWithoutScheme = s"$bucket/hello/world/with non ascii £€",
      parent = s"s3://$bucket/hello/world",
      getParent = s"s3://$bucket/hello/world",
      root = s"s3://$bucket/",
      name = "with non ascii £€",
      getFileName = s"s3://$bucket/with non ascii £€",
      getNameCount = 3,
      isAbsolute = true),

    // GoodPath(
    //   description = "a s3 uri path with encoded characters",
    //   path = s"s3://$bucket/hello/world/encoded                    paces",
    //   normalize = false,
    //   pathAsString = s"s3://$bucket/hello/world/encoded                    paces",
    //   pathWithoutScheme = s"$bucket/hello/world/encoded                    paces",
    //   parent = s"s3://$bucket/hello/world/",
    //   getParent = s"s3://$bucket/hello/world/",
    //   root = s"s3://$bucket/",
    //   name = "encoded                    paces",
    //   getFileName = s"s3://$bucket/encoded                    paces",
    //   getNameCount = 3,
    //   isAbsolute = true),

    GoodPath(
      description = "a s3 uri path with encoded characters",
      path = s"s3://$bucket/hello/world/encoded                    paces",
      normalize = false,
      pathAsString = s"s3://$bucket/hello/world/encoded                    paces",
      pathWithoutScheme = s"$bucket/hello/world/encoded                    paces",
      parent = s"s3://$bucket/hello/world",
      getParent = s"s3://$bucket/hello/world",
      root = s"s3://$bucket/",
      name = "encoded                    paces",
      getFileName = s"s3://$bucket/encoded                    paces",
      getNameCount = 3,
      isAbsolute = true),

    // TODO: In order for this to pass tests, S3Path needs to implement the
    //       Path trait directly and cannot inherit. We will work on this later
    // GoodPath(
    //   description = "a bucket only path",
    //   path = s"s3://$bucket",
    //   normalize = false,
    //   pathAsString = s"s3://$bucket/",
    //   pathWithoutScheme = s"$bucket/",
    //   parent = null,
    //   getParent = null,
    //   root = s"s3://$bucket/",
    //   name = "",
    //   getFileName = s"s3://$bucket/",
    //   getNameCount = 1,
    //   isAbsolute = false),

    GoodPath(
      description = "a bucket only path ending in a /",
      path = s"s3://$bucket/",
      normalize = false,
      pathAsString = s"s3://$bucket/",
      pathWithoutScheme = s"$bucket/",
      parent = null,
      getParent = null,
      root = s"s3://$bucket/",
      name = "",
      getFileName = null,
      getNameCount = 0,
      isAbsolute = true),

    GoodPath(
      description = "a file at the top of the bucket",
      path = s"s3://$bucket/hello",
      normalize = false,
      pathAsString = s"s3://$bucket/hello",
      pathWithoutScheme = s"$bucket/hello",
      parent = s"s3://$bucket/",
      getParent = s"s3://$bucket/",
      root = s"s3://$bucket/",
      name = "hello",
      getFileName = s"s3://$bucket/hello",
      getNameCount = 1,
      isAbsolute = true),

    // parent/getParent do not end in a "/".
    // TODO: Determine if this is critcal. Note
    //       that there are no "directories" in S3
    // GoodPath(
    //   description = "a path ending in /",
    //   path = s"s3://$bucket/hello/world/",
    //   normalize = false,
    //   pathAsString = s"s3://$bucket/hello/world/",
    //   pathWithoutScheme = s"$bucket/hello/world/",
    //   parent = s"s3://$bucket/hello/",
    //   getParent = s"s3://$bucket/hello/",
    //   root = s"s3://$bucket/",
    //   name = "world",
    //   getFileName = s"s3://$bucket/world",
    //   getNameCount = 2,
    //   isAbsolute = true),
    GoodPath(
      description = "a path ending in /",
      path = s"s3://$bucket/hello/world/",
      normalize = false,
      pathAsString = s"s3://$bucket/hello/world",
      pathWithoutScheme = s"$bucket/hello/world",
      parent = s"s3://$bucket/hello",
      getParent = s"s3://$bucket/hello",
      root = s"s3://$bucket/",
      name = "world",
      getFileName = s"s3://$bucket/world",
      getNameCount = 2,
      isAbsolute = true),

    // Special paths

    GoodPath(
      description = "a bucket with a path .",
      path = s"s3://$bucket/.",
      normalize = false,
      pathAsString = s"s3://$bucket/.",
      pathWithoutScheme = s"$bucket/.",
      parent = null,
      getParent = s"s3://$bucket/",
      root = s"s3://$bucket/",
      name = "",
      getFileName = s"s3://$bucket/.",
      getNameCount = 1,
      isAbsolute = true),

    // GoodPath(
    //   description = "a bucket with a path ..",
    //   path = s"s3://$bucket/..",
    //   normalize = false,
    //   pathAsString = s"s3://$bucket/..",
    //   pathWithoutScheme = s"$bucket/..",
    //   parent = null,
    //   getParent = s"s3://$bucket/",
    //   root = null,
    //   name = "",
    //   getFileName = s"s3://$bucket/..",
    //   getNameCount = 1,
    //   isAbsolute = true),

    // GoodPath(
    //   description = "a bucket including . in the path",
    //   path = s"s3://$bucket/hello/./world",
    //   normalize = false,
    //   pathAsString = s"s3://$bucket/hello/./world",
    //   pathWithoutScheme = s"$bucket/hello/./world",
    //   parent = s"s3://$bucket/hello/",
    //   getParent = s"s3://$bucket/hello/./",
    //   root = s"s3://$bucket/",
    //   name = "world",
    //   getFileName = s"s3://$bucket/world",
    //   getNameCount = 3,
    //   isAbsolute = true),
    //
    // GoodPath(
    //   description = "a bucket including .. in the path",
    //   path = s"s3://$bucket/hello/../world",
    //   normalize = false,
    //   pathAsString = s"s3://$bucket/hello/../world",
    //   pathWithoutScheme = s"$bucket/hello/../world",
    //   parent = s"s3://$bucket/",
    //   getParent = s"s3://$bucket/hello/../",
    //   root = s"s3://$bucket/",
    //   name = "world",
    //   getFileName = s"s3://$bucket/world",
    //   getNameCount = 3,
    //   isAbsolute = true),

    // Normalized

    // GoodPath(
    //   description = "a bucket with a normalized path .",
    //   path = s"s3://$bucket/.",
    //   normalize = true,
    //   pathAsString = s"s3://$bucket/",
    //   pathWithoutScheme = s"$bucket/",
    //   parent = null,
    //   getParent = null,
    //   root = s"s3://$bucket/",
    //   name = "",
    //   getFileName = null,
    //   getNameCount = 0,
    //   isAbsolute = true),
    //
    // GoodPath(
    //   description = "a bucket with a normalized path ..",
    //   path = s"s3://$bucket/..",
    //   normalize = true,
    //   pathAsString = s"s3://$bucket/",
    //   pathWithoutScheme = s"$bucket/",
    //   parent = null,
    //   getParent = null,
    //   root = s"s3://$bucket/",
    //   name = "",
    //   getFileName = s"s3://$bucket/",
    //   getNameCount = 1,
    //   isAbsolute = false),
    //
    // GoodPath(
    //   description = "a bucket including . in the normalized path",
    //   path = s"s3://$bucket/hello/./world",
    //   normalize = true,
    //   pathAsString = s"s3://$bucket/hello/world",
    //   pathWithoutScheme = s"$bucket/hello/world",
    //   parent = s"s3://$bucket/hello/",
    //   getParent = s"s3://$bucket/hello/",
    //   root = s"s3://$bucket/",
    //   name = "world",
    //   getFileName = s"s3://$bucket/world",
    //   getNameCount = 2,
    //   isAbsolute = true),

    GoodPath(
      description = "a bucket including .. in the normalized path",
      path = s"s3://$bucket/hello/../world",
      normalize = true,
      pathAsString = s"s3://$bucket/world",
      pathWithoutScheme = s"$bucket/world",
      parent = s"s3://$bucket/",
      getParent = s"s3://$bucket/",
      root = s"s3://$bucket/",
      name = "world",
      getFileName = s"s3://$bucket/world",
      getNameCount = 1,
      isAbsolute = true),

    GoodPath(
      description = "a bucket with an underscore",
      path = s"s3://hello_underscore/world",
      normalize = true,
      pathAsString = s"s3://hello_underscore/world",
      pathWithoutScheme = s"hello_underscore/world",
      parent = s"s3://hello_underscore/",
      getParent = s"s3://hello_underscore/",
      root = s"s3://hello_underscore/",
      name = "world",
      getFileName = s"s3://hello_underscore/world",
      getNameCount = 1,
      isAbsolute = true)
  )

  private def badPaths = Seq(
    BadPath("an empty path", "", " does not have a s3 scheme"),
    BadPath("a bucketless path", "s3://", "The specified S3 path 's3://' does not parse as a URI.\nExpected authority at index 5: s3://"),
    BadPath("a bucket named .", "s3://./hello/world", "The path 's3://./hello/world' does not seem to be a valid S3 path. Please check that it starts with s3:// and that the bucket and object follow S3 naming guidelines at https://docs.aws.amazon.com/AmazonS3/latest/dev/BucketRestrictions.html"),
    BadPath("a non ascii bucket name", "s3://nonasciibucket£€/hello/world",
      "The path 's3://nonasciibucket£€/hello/world' does not seem to be a valid S3 path. Please check that it starts with s3:// and that the bucket and object follow S3 naming guidelines at https://docs.aws.amazon.com/AmazonS3/latest/dev/BucketRestrictions.html"),
    BadPath("a https path", "https://hello/world", "S3 URIs must have 's3' scheme: https://hello/world"),
    BadPath("a file uri path", "file:///hello/world", "S3 URIs must have 's3' scheme: file:///hello/world"),
    BadPath("a relative file path", "hello/world", "hello/world does not have a s3 scheme"),
    BadPath("an absolute file path", "/hello/world", "/hello/world does not have a s3 scheme")
  )

  private lazy val pathBuilder = {
    S3PathBuilder.fromCredentials(
      AnonymousCredentialsProvider.create.resolveCredentials(),
      S3Storage.s3Configuration(),
      WorkflowOptions.empty,
      Option(Region.US_EAST_1)
    )
  }
}
