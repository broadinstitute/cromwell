package cromwell.core.path

import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpecLike, Matchers, Suite}

class MappedPathBuilderSpec extends Suite  with FlatSpecLike with Matchers with PathBuilderSpecUtils {
  behavior of "MappedPathBuilder"

  it should behave like truncateCommonRoots(prefixedPathBuilder, pathsToTruncate)

  goodPaths foreach { goodPath =>
    it should behave like buildGoodPath(prefixedPathBuilder, goodPath)
  }

  badPaths foreach { badPath =>
    it should behave like buildBadPath(prefixedPathBuilder, badPath)
  }

  private def pathsToTruncate = Table(
    ("context", "file", "relative"),
    ("gs://bucket", "gs://bucket/path/to/file", "path/to/file"),
    ("gs://bucket/path/to/my/dir", "gs://bucket/path/to/my/dir/file", "file"),
    ("gs://bucket/path/to/my/dir", "gs://bucket/path/to/my/dir//file", "file"),
    ("gs://bucket/path/to/my//dir", "gs://bucket/path/to/my/dir/file", "file"),
    ("gs://bucket/path/to/my//dir", "gs://bucket/path/to/my/dir//file", "file"),
    ("gs://bucket/path/to/my/dir", "gs://bucket/path/./to/my/dir/file", "./to/my/dir/file"),
    ("gs://bucket/path/to/my/dir/with/file", "gs://bucket/path/to/other/dir/with/file", "other/dir/with/file")
  )

  private def goodPaths = Seq(
    GoodPath(
      description = "a prefixed path",
      path = "gs://hello/world",
      normalize = false,
      pathAsString = "/my/mapped/dir/hello/world",
      pathWithoutScheme = "hello/world",
      parent = "/my/mapped/dir/hello",
      getParent = "/my/mapped/dir/hello",
      root = "/",
      name = "world",
      getFileName = "world",
      toUriHost = null,
      toUriPath = "/my/mapped/dir/hello/world",
      toUriStartsWith = "file:///",
      toUriEndsWith = "/my/mapped/dir/hello/world",
      getNameCount = 5,
      isAbsolute = true,
      isDirectory = false),

    GoodPath(
      description = "an empty prefixed path",
      path = "gs://",
      normalize = false,
      pathAsString = "/my/mapped/dir",
      pathWithoutScheme = "",
      parent = "/my/mapped",
      getParent = "/my/mapped",
      root = "/",
      name = "dir",
      getFileName = "dir",
      toUriHost = null,
      toUriPath = "/my/mapped/dir",
      toUriStartsWith = "file:///",
      toUriEndsWith = "/my/mapped/dir",
      getNameCount = 3,
      isAbsolute = true,
      isDirectory = false),

    GoodPath(
      description = "a prefixed path starting with .",
      path = "gs://./hello/world",
      normalize = false,
      pathAsString = "/my/mapped/dir/./hello/world",
      pathWithoutScheme = "./hello/world",
      parent = "/my/mapped/dir/hello",
      getParent = "/my/mapped/dir/./hello",
      root = "/",
      name = "world",
      getFileName = "world",
      toUriHost = null,
      toUriPath = "/my/mapped/dir/./hello/world",
      toUriStartsWith = "file:///",
      toUriEndsWith = "/my/mapped/dir/./hello/world",
      getNameCount = 6,
      isAbsolute = true,
      isDirectory = false)
  )

  private def badPaths = Seq(
    BadPath("an empty path", "", " must start with gs://"),
    BadPath("a https path", "https://hello/world", "https://hello/world must start with gs://"),
    BadPath("a file uri path", "file:///hello/world", "file:///hello/world must start with gs://"),
    BadPath("a relative file path", "hello/world", "hello/world must start with gs://"),
    BadPath("an absolute file path", "/hello/world", "/hello/world must start with gs://"),
    BadPath("an mapped directory path", "/my/mapped/dir/hello/world", "/my/mapped/dir/hello/world must start with gs://")
  )

  private lazy val prefixedPathBuilder = {
    new MappedPathBuilder("gs://", "/my/mapped/dir")
  }
}
