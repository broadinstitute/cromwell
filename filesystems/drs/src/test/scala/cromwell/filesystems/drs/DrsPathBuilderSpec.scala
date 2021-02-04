package cromwell.filesystems.drs

import cloud.nio.impl.drs.DrsCloudNioFileProvider.DrsReadInterpreter
import cloud.nio.impl.drs.DrsCloudNioFileSystemProvider
import com.google.cloud.NoCredentials
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.TestKitSuite
import cromwell.core.path._
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.Tables.Table

class DrsPathBuilderSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with PathBuilderSpecUtils {

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
      parent = null,
      getParent = null,
      root = s"drs://$bucket/hello/world/with spaces",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    GoodPath(
      description = "a path with non-ascii",
      path = s"drs://$bucket/hello/world/with non ascii £€",
      normalize = false,
      pathAsString = s"drs://$bucket/hello/world/with non ascii £€",
      pathWithoutScheme = s"$bucket/hello/world/with non ascii £€",
      parent = null,
      getParent = null,
      root = s"drs://$bucket/hello/world/with non ascii £€",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    GoodPath(
      description = "a gs uri path with encoded characters",
      path = s"drs://$bucket/hello/world/encoded%20spaces",
      normalize = false,
      pathAsString = s"drs://$bucket/hello/world/encoded%20spaces",
      pathWithoutScheme = s"$bucket/hello/world/encoded%20spaces",
      parent = null,
      getParent = null,
      root = s"drs://$bucket/hello/world/encoded%20spaces",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    GoodPath(
      description = "a file at the top of the bucket",
      path = s"drs://$bucket/hello",
      normalize = false,
      pathAsString = s"drs://$bucket/hello",
      pathWithoutScheme = s"$bucket/hello",
      parent = null,
      getParent = null,
      root = s"drs://$bucket/hello",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    GoodPath(
      description = "a path ending in /",
      path = s"drs://$bucket/hello/world/",
      normalize = false,
      pathAsString = s"drs://$bucket/hello/world/",
      pathWithoutScheme = s"$bucket/hello/world/",
      parent = null,
      getParent = null,
      root = s"drs://$bucket/hello/world/",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    // Special paths

    GoodPath(
      description = "a bucket with a path .",
      path = s"drs://$bucket/.",
      normalize = false,
      pathAsString = s"drs://$bucket/.",
      pathWithoutScheme = s"$bucket/.",
      parent = null,
      getParent = null,
      root = s"drs://$bucket/.",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    GoodPath(
      description = "a bucket with a path ..",
      path = s"drs://$bucket/..",
      normalize = false,
      pathAsString = s"drs://$bucket/..",
      pathWithoutScheme = s"$bucket/..",
      parent = null,
      getParent = null,
      root = s"drs://$bucket/..",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    GoodPath(
      description = "a bucket including . in the path",
      path = s"drs://$bucket/hello/./world",
      normalize = false,
      pathAsString = s"drs://$bucket/hello/./world",
      pathWithoutScheme = s"$bucket/hello/./world",
      parent = null,
      getParent = null,
      root = s"drs://$bucket/hello/./world",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    GoodPath(
      description = "a bucket including .. in the path",
      path = s"drs://$bucket/hello/../world",
      normalize = false,
      pathAsString = s"drs://$bucket/hello/../world",
      pathWithoutScheme = s"$bucket/hello/../world",
      parent = null,
      getParent = null,
      root = s"drs://$bucket/hello/../world",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    // Normalized

    GoodPath(
      description = "a bucket with a normalized path .",
      path = s"drs://$bucket/.",
      normalize = true,
      pathAsString = s"drs://$bucket/.",
      pathWithoutScheme = s"$bucket/.",
      parent = null,
      getParent = null,
      root = s"drs://$bucket/.",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    GoodPath(
      description = "a bucket with a normalized path ..",
      path = s"drs://$bucket/..",
      normalize = true,
      pathAsString = s"drs://$bucket/..",
      pathWithoutScheme = s"$bucket/..",
      parent = null,
      getParent = null,
      root = s"drs://$bucket/..",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    GoodPath(
      description = "a bucket including . in the normalized path",
      path = s"drs://$bucket/hello/./world",
      normalize = true,
      pathAsString = s"drs://$bucket/hello/./world",
      pathWithoutScheme = s"$bucket/hello/./world",
      parent = null,
      getParent = null,
      root = s"drs://$bucket/hello/./world",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    GoodPath(
      description = "a bucket including .. in the normalized path",
      path = s"drs://$bucket/hello/../world",
      normalize = true,
      pathAsString = s"drs://$bucket/hello/../world",
      pathWithoutScheme = s"$bucket/hello/../world",
      parent = null,
      getParent = null,
      root = s"drs://$bucket/hello/../world",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    GoodPath(
      description = "a bucket with an underscore",
      path = s"drs://hello_underscore/world",
      normalize = true,
      pathAsString = s"drs://hello_underscore/world",
      pathWithoutScheme = s"hello_underscore/world",
      parent = null,
      getParent = null,
      root = s"drs://hello_underscore/world",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    GoodPath(
      description = "a bucket named .",
      path = s"drs://./hello/world",
      normalize = true,
      pathAsString = s"drs://./hello/world",
      pathWithoutScheme = s"./hello/world",
      parent = null,
      getParent = null,
      root = s"drs://./hello/world",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    GoodPath(
      description = "a non ascii bucket name",
      path = s"drs://nonasciibucket£€/hello/world",
      normalize = true,
      pathAsString = s"drs://nonasciibucket£€/hello/world",
      pathWithoutScheme = s"nonasciibucket£€/hello/world",
      parent = null,
      getParent = null,
      root = s"drs://nonasciibucket£€/hello/world",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

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
      isAbsolute = false,
    ),

    GoodPath(
      description = "an absolute path without a host",
      path = s"drs://blah",
      normalize = false,
      pathAsString = s"drs://blah",
      pathWithoutScheme = s"blah",
      parent = null,
      getParent = null,
      root = s"drs://blah",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    // No spec says this is illegal... so pass it to Martha's various GCFs JIC
    GoodPath(
      description = "a bucketless path",
      path = s"drs://",
      normalize = false,
      pathAsString = s"drs://",
      pathWithoutScheme = "",
      parent = null,
      getParent = null,
      root = s"drs://",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    // Sample via: https://docs.google.com/document/d/1Wf4enSGOEXD5_AE-uzLoYqjIp5MnePbZ6kYTVFp1WoM/edit
    GoodPath(
      description = "a path with a query string",
      path = "drs://drs.data.humancellatlas.org/8aca942c-17f7-4e34-b8fd-3c12e50f9291?version=2019-07-04T151444.185805Z",
      normalize = false,
      pathAsString =
        "drs://drs.data.humancellatlas.org/8aca942c-17f7-4e34-b8fd-3c12e50f9291?version=2019-07-04T151444.185805Z",
      pathWithoutScheme =
        "drs.data.humancellatlas.org/8aca942c-17f7-4e34-b8fd-3c12e50f9291?version=2019-07-04T151444.185805Z",
      parent = null,
      getParent = null,
      root =
        "drs://drs.data.humancellatlas.org/8aca942c-17f7-4e34-b8fd-3c12e50f9291?version=2019-07-04T151444.185805Z",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),

    // Sample via: https://docs.google.com/document/d/1Wf4enSGOEXD5_AE-uzLoYqjIp5MnePbZ6kYTVFp1WoM/edit
    GoodPath(
      description = "a compact identifier-based path",
      path = s"drs://dg.ANV0:dg.ANV0/0db6577e-57bd-48a1-93c6-327c292bcb6b",
      normalize = false,
      pathAsString = s"drs://dg.ANV0:dg.ANV0/0db6577e-57bd-48a1-93c6-327c292bcb6b",
      pathWithoutScheme = "dg.ANV0:dg.ANV0/0db6577e-57bd-48a1-93c6-327c292bcb6b",
      parent = null,
      getParent = null,
      root = s"drs://dg.ANV0:dg.ANV0/0db6577e-57bd-48a1-93c6-327c292bcb6b",
      name = "",
      getFileName = null,
      getNameCount = 1,
      isAbsolute = false,
    ),
  )

  private def badPaths = Seq(
    BadPath("an empty path", "", " does not have a drs scheme."),
    BadPath("a GCS path", s"gs://$bucket/hello/world", "gs://mymadeupbucket/hello/world does not have a drs scheme."),
    BadPath("a https path", "https://hello/world", "https://hello/world does not have a drs scheme."),
    BadPath("a file uri path", "file:///hello/world", "file:///hello/world does not have a drs scheme."),
    BadPath("a relative file path", "hello/world", "hello/world does not have a drs scheme."),
    BadPath("an absolute file path", "/hello/world", "/hello/world does not have a drs scheme."),
  )

  private val drsReadInterpreter: DrsReadInterpreter = (_, _) =>
    throw new UnsupportedOperationException("Currently DrsPathBuilderSpec doesn't need to use drs read interpreter.")

  private val marthaConfig: Config = ConfigFactory.parseString(
    """martha {
      |   url = "http://martha-url"
      |}
      |""".stripMargin
  )

  private lazy val fakeCredentials = NoCredentials.getInstance

  private lazy val drsPathBuilder = DrsPathBuilder(
    new DrsCloudNioFileSystemProvider(marthaConfig, fakeCredentials, drsReadInterpreter),
    None,
  )
}
