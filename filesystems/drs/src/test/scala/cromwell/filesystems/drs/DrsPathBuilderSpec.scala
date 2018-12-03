package cromwell.filesystems.drs

import cloud.nio.impl.drs.DrsCloudNioFileSystemProvider
import com.typesafe.config.ConfigFactory
import cromwell.core.TestKitSuite
import cromwell.core.path._
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
      ("dos://bucket", "dos://bucket/path/to/file", "path/to/file"),
      ("dos://bucket/path/to/my/dir", "dos://bucket/path/to/my/dir/file", "file"),
      ("dos://bucket/path/to/my/dir", "dos://bucket/path/to/my/dir//file", "file"),
      // NOTE: Next two are different from the DefaultPathBuilder. "//" doesn't build to "/" in the GcsPathBuilder
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
        description = "a bucket only path",
        path = s"dos://$bucket",
        normalize = false,
        pathAsString = s"dos://$bucket/",
        pathWithoutScheme = s"$bucket/",
        parent = null,
        getParent = null,
        root = s"dos://$bucket/",
        name = "",
        getFileName = s"dos://$bucket/",
        getNameCount = 1,
        isAbsolute = false),

      GoodPath(
        description = "a bucket only path ending in a /",
        path = s"dos://$bucket/",
        normalize = false,
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
        pathAsString = s"dos://$bucket/hello/world/",
        pathWithoutScheme = s"$bucket/hello/world/",
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
        /*
        Begin DOS/GCS differences after normalize

        Has to do with how the paths are converted to strings and then re-parsed. After re-parsing the internal `"../"` is
        lost and converted to `""`.
        */
        getFileName = null,
        getNameCount = 0,
        isAbsolute = true),
        /*
        End DOS/GCS differences after normalize
        */

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
        isAbsolute = true)
    )

    private def badPaths = Seq(
      BadPath("an empty path", "", " does not have a dos scheme"),
      BadPath("a GCS path", s"gs://$bucket/hello/world", "gs://mymadeupbucket/hello/world does not have a dos scheme"),
      BadPath("an bucketless path", "dos://", "The specified DOS path 'dos://' does not parse as a URI.\nExpected authority at index 5: dos://"),
      BadPath("a bucket named .", "dos://./hello/world", "The path 'dos://./hello/world' does not seem to be a valid DOS path. Please check that it starts with dos:// and that the bucket and object follow DOS naming guidelines at https://cloud.google.com/storage/docs/naming."),
      BadPath("a non ascii bucket name", "dos://nonasciibucket£€/hello/world",
        "The path 'dos://nonasciibucket£€/hello/world' does not seem to be a valid DOS path. Please check that it starts with dos:// and that the bucket and object follow DOS naming guidelines at https://cloud.google.com/storage/docs/naming."),
      BadPath("a https path", "https://hello/world", "https://hello/world does not have a dos scheme"),
      BadPath("a file uri path", "file:///hello/world", "file:///hello/world does not have a dos scheme"),
      BadPath("a relative file path", "hello/world", "hello/world does not have a dos scheme"),
      BadPath("an absolute file path", "/hello/world", "/hello/world does not have a dos scheme")
    )

  val globalConfig = ConfigFactory.parseString(
    """|filesystems {
       |  drs {
       |    class = "cromwell.filesystems.drs.DrsPathBuilderFactory"
       |    global {
       |      class = "cromwell.filesystems.drs.DrsFileSystemConfig"
       |      config {
       |        martha {
       |          url = "http://matha-url"
       |          request.json-template = "{"key": "${holder}"}"
       |        }
       |      }
       |    }
       |  }
       |}
       |""".stripMargin
  )

    private lazy val drsPathBuilder = DrsPathBuilder(new DrsCloudNioFileSystemProvider(globalConfig))
}
