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
    ("drs://bucket/path/to/my/dir", "drs://bucket/path/to/my/dir/file", "file"),
    ("drs://bucket/path/to/my/dir", "drs://bucket/path/to/my/dir//file", "file"),
    ("drs://bucket/path/to/my//dir", "drs://bucket/path/to/my/dir/file", "dir/file"),
    ("drs://bucket/path/to/my//dir", "drs://bucket/path/to/my/dir//file", "dir//file"),
    ("drs://bucket/path/to/my/dir", "drs://bucket/path/./to/my/dir/file", "./to/my/dir/file"),
    ("drs://bucket/path/to/my/dir/with/file", "drs://bucket/path/to/other/dir/with/file", "other/dir/with/file")
  )

  private def bucket = "mymadeupbucket"

  private def goodPaths = Seq(
    GoodPath(
      description = "a path with spaces",
      path = s"drs://$bucket/hello/world/with spaces",
      normalize = false,
      pathAsString = s"drs://$bucket/hello/world/with spaces",
      pathWithoutScheme = s"$bucket/hello/world/with spaces",
      parent = s"drs://$bucket/hello/world/",
      getParent = s"drs://$bucket/hello/world/",
      root = s"drs://$bucket/",
      name = "with spaces",
      getFileName = s"drs://$bucket/with spaces",
      getNameCount = 3,
      isAbsolute = true),

    GoodPath(
      description = "a path with non-ascii",
      path = s"drs://$bucket/hello/world/with non ascii £€",
      normalize = false,
      pathAsString = s"drs://$bucket/hello/world/with non ascii £€",
      pathWithoutScheme = s"$bucket/hello/world/with non ascii £€",
      parent = s"drs://$bucket/hello/world/",
      getParent = s"drs://$bucket/hello/world/",
      root = s"drs://$bucket/",
      name = "with non ascii £€",
      getFileName = s"drs://$bucket/with non ascii £€",
      getNameCount = 3,
      isAbsolute = true),

    GoodPath(
      description = "a gs uri path with encoded characters",
      path = s"drs://$bucket/hello/world/encoded%20spaces",
      normalize = false,
      pathAsString = s"drs://$bucket/hello/world/encoded%20spaces",
      pathWithoutScheme = s"$bucket/hello/world/encoded%20spaces",
      parent = s"drs://$bucket/hello/world/",
      getParent = s"drs://$bucket/hello/world/",
      root = s"drs://$bucket/",
      name = "encoded%20spaces",
      getFileName = s"drs://$bucket/encoded%20spaces",
      getNameCount = 3,
      isAbsolute = true),

    GoodPath(
      description = "a file at the top of the bucket",
      path = s"drs://$bucket/hello",
      normalize = false,
      pathAsString = s"drs://$bucket/hello",
      pathWithoutScheme = s"$bucket/hello",
      parent = s"drs://$bucket/",
      getParent = s"drs://$bucket/",
      root = s"drs://$bucket/",
      name = "hello",
      getFileName = s"drs://$bucket/hello",
      getNameCount = 1,
      isAbsolute = true),

    GoodPath(
      description = "a path ending in /",
      path = s"drs://$bucket/hello/world/",
      normalize = false,
      pathAsString = s"drs://$bucket/hello/world",
      pathWithoutScheme = s"$bucket/hello/world",
      parent = s"drs://$bucket/hello/",
      getParent = s"drs://$bucket/hello/",
      root = s"drs://$bucket/",
      name = "world",
      getFileName = s"drs://$bucket/world",
      getNameCount = 2,
      isAbsolute = true),

    // Special paths

    GoodPath(
      description = "a bucket with a path .",
      path = s"drs://$bucket/.",
      normalize = false,
      pathAsString = s"drs://$bucket/.",
      pathWithoutScheme = s"$bucket/.",
      parent = null,
      getParent = s"drs://$bucket/",
      root = s"drs://$bucket/",
      name = "",
      getFileName = s"drs://$bucket/.",
      getNameCount = 1,
      isAbsolute = true),

    GoodPath(
      description = "a bucket with a path ..",
      path = s"drs://$bucket/..",
      normalize = false,
      pathAsString = s"drs://$bucket/..",
      pathWithoutScheme = s"$bucket/..",
      parent = null,
      getParent = s"drs://$bucket/",
      root = null,
      name = "",
      getFileName = s"drs://$bucket/..",
      getNameCount = 1,
      isAbsolute = true),

    GoodPath(
      description = "a bucket including . in the path",
      path = s"drs://$bucket/hello/./world",
      normalize = false,
      pathAsString = s"drs://$bucket/hello/./world",
      pathWithoutScheme = s"$bucket/hello/./world",
      parent = s"drs://$bucket/hello/",
      getParent = s"drs://$bucket/hello/./",
      root = s"drs://$bucket/",
      name = "world",
      getFileName = s"drs://$bucket/world",
      getNameCount = 3,
      isAbsolute = true),

    GoodPath(
      description = "a bucket including .. in the path",
      path = s"drs://$bucket/hello/../world",
      normalize = false,
      pathAsString = s"drs://$bucket/hello/../world",
      pathWithoutScheme = s"$bucket/hello/../world",
      parent = s"drs://$bucket/",
      getParent = s"drs://$bucket/hello/../",
      root = s"drs://$bucket/",
      name = "world",
      getFileName = s"drs://$bucket/world",
      getNameCount = 3,
      isAbsolute = true),

    // Normalized

    GoodPath(
      description = "a bucket with a normalized path .",
      path = s"drs://$bucket/.",
      normalize = true,
      pathAsString = s"drs://$bucket/",
      pathWithoutScheme = s"$bucket/",
      parent = null,
      getParent = null,
      root = s"drs://$bucket/",
      name = "",
      getFileName = null,
      getNameCount = 0,
      isAbsolute = true),

    GoodPath(
      description = "a bucket with a normalized path ..",
      path = s"drs://$bucket/..",
      normalize = true,
      pathAsString = s"drs://$bucket/",
      pathWithoutScheme = s"$bucket/",
      parent = null,
      getParent = null,
      root = s"drs://$bucket/",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false),

    GoodPath(
      description = "a bucket including . in the normalized path",
      path = s"drs://$bucket/hello/./world",
      normalize = true,
      pathAsString = s"drs://$bucket/hello/world",
      pathWithoutScheme = s"$bucket/hello/world",
      parent = s"drs://$bucket/hello/",
      getParent = s"drs://$bucket/hello/",
      root = s"drs://$bucket/",
      name = "world",
      getFileName = s"drs://$bucket/world",
      getNameCount = 2,
      isAbsolute = true),

    GoodPath(
      description = "a bucket including .. in the normalized path",
      path = s"drs://$bucket/hello/../world",
      normalize = true,
      pathAsString = s"drs://$bucket/world",
      pathWithoutScheme = s"$bucket/world",
      parent = s"drs://$bucket/",
      getParent = s"drs://$bucket/",
      root = s"drs://$bucket/",
      name = "world",
      getFileName = s"drs://$bucket/world",
      getNameCount = 1,
      isAbsolute = true),

    GoodPath(
      description = "a bucket with an underscore",
      path = s"drs://hello_underscore/world",
      normalize = true,
      pathAsString = s"drs://hello_underscore/world",
      pathWithoutScheme = s"hello_underscore/world",
      parent = s"drs://hello_underscore/",
      getParent = s"drs://hello_underscore/",
      root = s"drs://hello_underscore/",
      name = "world",
      getFileName = s"drs://hello_underscore/world",
      getNameCount = 1,
      isAbsolute = true),

    GoodPath(
      description = "a bucket named .",
      path = s"drs://./hello/world",
      normalize = true,
      pathAsString = s"drs://./hello/world",
      pathWithoutScheme = s"./hello/world",
      parent = s"drs://./hello/",
      getParent = s"drs://./hello/",
      root = s"drs://./",
      name = "world",
      getFileName = s"drs://./world",
      getNameCount = 2,
      isAbsolute = true),

    GoodPath(
      description = "a non ascii bucket name",
      path = s"drs://nonasciibucket£€/hello/world",
      normalize = true,
      pathAsString = s"drs://nonasciibucket£€/hello/world",
      pathWithoutScheme = s"nonasciibucket£€/hello/world",
      parent = s"drs://nonasciibucket£€/hello/",
      getParent = s"drs://nonasciibucket£€/hello/",
      root = s"drs://nonasciibucket£€/",
      name = "world",
      getFileName = s"drs://nonasciibucket£€/world",
      getNameCount = 2,
      isAbsolute = true),

    GoodPath(
      description = "an non-absolute path without a host",
      path = s"drs://blah/",
      normalize = false,
      pathAsString = s"drs://blah/",
      pathWithoutScheme = s"blah/",
      parent = null,
      getParent = null,
      root = s"drs://blah/",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false),

    GoodPath(
      description = "an absolute path without a host",
      path = s"drs://blah",
      normalize = false,
      pathAsString = s"drs://blah/",
      pathWithoutScheme = s"blah/",
      parent = null,
      getParent = null,
      root = s"drs://blah/",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false)
  )

  private def badPaths = Seq(
    BadPath("an empty path", "", " does not have a drs scheme."),
    BadPath("a GCS path", s"gs://$bucket/hello/world", "gs://mymadeupbucket/hello/world does not have a drs scheme."),
    BadPath("an bucketless path", "drs://", "Expected authority at index 6: drs://"),
    BadPath("a https path", "https://hello/world", "https://hello/world does not have a drs scheme."),
    BadPath("a file uri path", "file:///hello/world", "file:///hello/world does not have a drs scheme."),
    BadPath("a relative file path", "hello/world", "hello/world does not have a drs scheme."),
    BadPath("an absolute file path", "/hello/world", "/hello/world does not have a drs scheme."),
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
