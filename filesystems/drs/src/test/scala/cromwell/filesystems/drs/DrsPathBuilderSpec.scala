package cromwell.filesystems.drs

import java.nio.channels.ReadableByteChannel

import cats.effect.IO
import cloud.nio.impl.drs.{DrsCloudNioFileSystemProvider, MarthaResponse}
import com.google.cloud.NoCredentials
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.TestKitSuite
import cromwell.core.path._
import org.apache.http.impl.client.HttpClientBuilder
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpecLike, Matchers}

class DrsPathBuilderSpec extends TestKitSuite with FlatSpecLike with Matchers with PathBuilderSpecUtils {

  behavior of "DrsPathBuilder"

  it should behave like truncateCommonRoots(drsPathBuilder, pathsToTruncate)

  goodPaths foreach { goodPath =>
    it should behave like buildGoodPath(drsPathBuilder, goodPath)
  }

  badPaths foreach { badPath =>
    it should behave like buildBadPath(drsPathBuilder, badPath)
  }

  private def pathsToTruncate = Table(
    ("context", "file", "relative"),
    ("dos://bucket/path/to/my/dir", "dos://bucket/path/to/my/dir/file", "file"),
    ("dos://bucket/path/to/my/dir", "dos://bucket/path/to/my/dir//file", "file"),
    ("dos://bucket/path/to/my//dir", "dos://bucket/path/to/my/dir/file", "dir/file"),
    ("dos://bucket/path/to/my//dir", "dos://bucket/path/to/my/dir//file", "dir//file"),
    ("dos://bucket/path/to/my/dir", "dos://bucket/path/./to/my/dir/file", "./to/my/dir/file"),
    ("dos://bucket/path/to/my/dir/with/file", "dos://bucket/path/to/other/dir/with/file", "other/dir/with/file")
  )

  private def bucket = "mymadeupbucket"

  private def goodPaths = Seq(
    GoodPath(
      description = "a path with spaces",
      path = s"dos://$bucket/hello/world/with spaces",
      normalize = false,
      pathAsString = s"dos://$bucket/hello/world/with spaces",
      pathWithoutScheme = s"$bucket/hello/world/with spaces",
      parent = s"dos://$bucket/hello/world/",
      getParent = s"dos://$bucket/hello/world/",
      root = s"dos://$bucket/",
      name = "with spaces",
      getFileName = s"dos://$bucket/with spaces",
      getNameCount = 3,
      isAbsolute = true),

    GoodPath(
      description = "a path with non-ascii",
      path = s"dos://$bucket/hello/world/with non ascii £€",
      normalize = false,
      pathAsString = s"dos://$bucket/hello/world/with non ascii £€",
      pathWithoutScheme = s"$bucket/hello/world/with non ascii £€",
      parent = s"dos://$bucket/hello/world/",
      getParent = s"dos://$bucket/hello/world/",
      root = s"dos://$bucket/",
      name = "with non ascii £€",
      getFileName = s"dos://$bucket/with non ascii £€",
      getNameCount = 3,
      isAbsolute = true),

    GoodPath(
      description = "a gs uri path with encoded characters",
      path = s"dos://$bucket/hello/world/encoded%20spaces",
      normalize = false,
      pathAsString = s"dos://$bucket/hello/world/encoded%20spaces",
      pathWithoutScheme = s"$bucket/hello/world/encoded%20spaces",
      parent = s"dos://$bucket/hello/world/",
      getParent = s"dos://$bucket/hello/world/",
      root = s"dos://$bucket/",
      name = "encoded%20spaces",
      getFileName = s"dos://$bucket/encoded%20spaces",
      getNameCount = 3,
      isAbsolute = true),

    GoodPath(
      description = "a file at the top of the bucket",
      path = s"dos://$bucket/hello",
      normalize = false,
      pathAsString = s"dos://$bucket/hello",
      pathWithoutScheme = s"$bucket/hello",
      parent = s"dos://$bucket/",
      getParent = s"dos://$bucket/",
      root = s"dos://$bucket/",
      name = "hello",
      getFileName = s"dos://$bucket/hello",
      getNameCount = 1,
      isAbsolute = true),

    GoodPath(
      description = "a path ending in /",
      path = s"dos://$bucket/hello/world/",
      normalize = false,
      pathAsString = s"dos://$bucket/hello/world",
      pathWithoutScheme = s"$bucket/hello/world",
      parent = s"dos://$bucket/hello/",
      getParent = s"dos://$bucket/hello/",
      root = s"dos://$bucket/",
      name = "world",
      getFileName = s"dos://$bucket/world",
      getNameCount = 2,
      isAbsolute = true),

    // Special paths

    GoodPath(
      description = "a bucket with a path .",
      path = s"dos://$bucket/.",
      normalize = false,
      pathAsString = s"dos://$bucket/.",
      pathWithoutScheme = s"$bucket/.",
      parent = null,
      getParent = s"dos://$bucket/",
      root = s"dos://$bucket/",
      name = "",
      getFileName = s"dos://$bucket/.",
      getNameCount = 1,
      isAbsolute = true),

    GoodPath(
      description = "a bucket with a path ..",
      path = s"dos://$bucket/..",
      normalize = false,
      pathAsString = s"dos://$bucket/..",
      pathWithoutScheme = s"$bucket/..",
      parent = null,
      getParent = s"dos://$bucket/",
      root = null,
      name = "",
      getFileName = s"dos://$bucket/..",
      getNameCount = 1,
      isAbsolute = true),

    GoodPath(
      description = "a bucket including . in the path",
      path = s"dos://$bucket/hello/./world",
      normalize = false,
      pathAsString = s"dos://$bucket/hello/./world",
      pathWithoutScheme = s"$bucket/hello/./world",
      parent = s"dos://$bucket/hello/",
      getParent = s"dos://$bucket/hello/./",
      root = s"dos://$bucket/",
      name = "world",
      getFileName = s"dos://$bucket/world",
      getNameCount = 3,
      isAbsolute = true),

    GoodPath(
      description = "a bucket including .. in the path",
      path = s"dos://$bucket/hello/../world",
      normalize = false,
      pathAsString = s"dos://$bucket/hello/../world",
      pathWithoutScheme = s"$bucket/hello/../world",
      parent = s"dos://$bucket/",
      getParent = s"dos://$bucket/hello/../",
      root = s"dos://$bucket/",
      name = "world",
      getFileName = s"dos://$bucket/world",
      getNameCount = 3,
      isAbsolute = true),

    // Normalized

    GoodPath(
      description = "a bucket with a normalized path .",
      path = s"dos://$bucket/.",
      normalize = true,
      pathAsString = s"dos://$bucket/",
      pathWithoutScheme = s"$bucket/",
      parent = null,
      getParent = null,
      root = s"dos://$bucket/",
      name = "",
      getFileName = null,
      getNameCount = 0,
      isAbsolute = true),

    GoodPath(
      description = "a bucket with a normalized path ..",
      path = s"dos://$bucket/..",
      normalize = true,
      pathAsString = s"dos://$bucket/",
      pathWithoutScheme = s"$bucket/",
      parent = null,
      getParent = null,
      root = s"dos://$bucket/",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false),

    GoodPath(
      description = "a bucket including . in the normalized path",
      path = s"dos://$bucket/hello/./world",
      normalize = true,
      pathAsString = s"dos://$bucket/hello/world",
      pathWithoutScheme = s"$bucket/hello/world",
      parent = s"dos://$bucket/hello/",
      getParent = s"dos://$bucket/hello/",
      root = s"dos://$bucket/",
      name = "world",
      getFileName = s"dos://$bucket/world",
      getNameCount = 2,
      isAbsolute = true),

    GoodPath(
      description = "a bucket including .. in the normalized path",
      path = s"dos://$bucket/hello/../world",
      normalize = true,
      pathAsString = s"dos://$bucket/world",
      pathWithoutScheme = s"$bucket/world",
      parent = s"dos://$bucket/",
      getParent = s"dos://$bucket/",
      root = s"dos://$bucket/",
      name = "world",
      getFileName = s"dos://$bucket/world",
      getNameCount = 1,
      isAbsolute = true),

    GoodPath(
      description = "a bucket with an underscore",
      path = s"dos://hello_underscore/world",
      normalize = true,
      pathAsString = s"dos://hello_underscore/world",
      pathWithoutScheme = s"hello_underscore/world",
      parent = s"dos://hello_underscore/",
      getParent = s"dos://hello_underscore/",
      root = s"dos://hello_underscore/",
      name = "world",
      getFileName = s"dos://hello_underscore/world",
      getNameCount = 1,
      isAbsolute = true),

    GoodPath(
      description = "a bucket named .",
      path = s"dos://./hello/world",
      normalize = true,
      pathAsString = s"dos://./hello/world",
      pathWithoutScheme = s"./hello/world",
      parent = s"dos://./hello/",
      getParent = s"dos://./hello/",
      root = s"dos://./",
      name = "world",
      getFileName = s"dos://./world",
      getNameCount = 2,
      isAbsolute = true),

    GoodPath(
      description = "a non ascii bucket name",
      path = s"dos://nonasciibucket£€/hello/world",
      normalize = true,
      pathAsString = s"dos://nonasciibucket£€/hello/world",
      pathWithoutScheme = s"nonasciibucket£€/hello/world",
      parent = s"dos://nonasciibucket£€/hello/",
      getParent = s"dos://nonasciibucket£€/hello/",
      root = s"dos://nonasciibucket£€/",
      name = "world",
      getFileName = s"dos://nonasciibucket£€/world",
      getNameCount = 2,
      isAbsolute = true),

    GoodPath(
      description = "an non-absolute path without a host",
      path = s"dos://blah/",
      normalize = false,
      pathAsString = s"dos://blah/",
      pathWithoutScheme = s"blah/",
      parent = null,
      getParent = null,
      root = s"dos://blah/",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false),

    GoodPath(
      description = "an absolute path without a host",
      path = s"dos://blah",
      normalize = false,
      pathAsString = s"dos://blah/",
      pathWithoutScheme = s"blah/",
      parent = null,
      getParent = null,
      root = s"dos://blah/",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false)
  )

  private def badPaths = Seq(
    BadPath("an empty path", "", " does not have a dos scheme."),
    BadPath("a GCS path", s"gs://$bucket/hello/world", "gs://mymadeupbucket/hello/world does not have a dos scheme."),
    BadPath("an bucketless path", "dos://", "Expected authority at index 6: dos://"),
    BadPath("a https path", "https://hello/world", "https://hello/world does not have a dos scheme."),
    BadPath("a file uri path", "file:///hello/world", "file:///hello/world does not have a dos scheme."),
    BadPath("a relative file path", "hello/world", "hello/world does not have a dos scheme."),
    BadPath("an absolute file path", "/hello/world", "/hello/world does not have a dos scheme."),
  )

  private def drsReadInterpreter(marthaResponse: MarthaResponse): IO[ReadableByteChannel] =
    throw new UnsupportedOperationException("Currently DrsPathBuilderSpec doesn't need to use drs read interpreter.")


  private val marthaConfig: Config = ConfigFactory.parseString(
    """martha {
      |   url = "http://martha-url"
      |   request.json-template = "{"key": "${holder}"}"
      |}
      |""".stripMargin
  )

  private lazy val fakeCredentials = NoCredentials.getInstance

  private lazy val httpClientBuilder = HttpClientBuilder.create()

  private lazy val drsPathBuilder = DrsPathBuilder(
    new DrsCloudNioFileSystemProvider(marthaConfig, fakeCredentials, httpClientBuilder, drsReadInterpreter),
    None,
  )
}
