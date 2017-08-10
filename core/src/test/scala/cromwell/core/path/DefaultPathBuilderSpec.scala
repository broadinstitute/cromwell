package cromwell.core.path

import cromwell.util.TestFileUtil
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpecLike, Matchers, Suite}

class DefaultPathBuilderSpec extends Suite with FlatSpecLike with Matchers with PathBuilderSpecUtils with TestFileUtil {

  private val pwd = BetterFileMethods.Cmds.pwd
  private val parentOption = Option(pwd.parent)
  private val grandParentOption = parentOption.flatMap(parent => Option(parent.parent))

  private val pwdPath = pwd.pathAsString
  private val pwdPathPrefix = if (pwd.pathAsString == "/") "" else pwdPath

  private val parentPath = parentOption.map(_.pathAsString).orNull
  private val grandParentPath = grandParentOption.map(_.pathAsString).orNull

  private val pwdFileNameOption = Option(pwd.getFileName)
  private val pwdName = pwdFileNameOption.map(_.pathAsString).getOrElse("")

  private val parentFileNameOption = parentOption.flatMap(parent => Option(parent.getFileName))
  private val parentName = parentFileNameOption.map(_.pathAsString).getOrElse("")

  behavior of "DefaultPathBuilder"

  it should "generate a md5 hash from a file" in {
    val f1 = createCannedFile("file1", "content")
    val f2 = createCannedFile("file2", "content")
    val f3 = createCannedFile("file3", "other content")
    f1.md5 should be("9A0364B9E99BB480DD25E1F0284C8555")
    f1.md5 should be(f2.md5)
    f1.md5 shouldNot be(f3.md5)
  }

  it should behave like truncateCommonRoots(DefaultPathBuilder, pathsToTruncate)

  goodPaths foreach { goodPath =>
    it should behave like buildGoodPath(DefaultPathBuilder, goodPath)
  }

  badPaths foreach { badPath =>
    it should behave like buildBadPath(DefaultPathBuilder, badPath)
  }

  private def pathsToTruncate = Table(
    ("context", "file", "relative"),
    ("/", "/path/to/file", "path/to/file"),
    ("/path/to/my/dir", "/path/to/my/dir/file", "file"),
    ("/path/to/my/dir", "/path/to/my/dir//file", "file"),
    ("/path/to/my//dir", "/path/to/my/dir/file", "file"),
    ("/path/to/my//dir", "/path/to/my/dir//file", "file"),
    ("/path/to/my/dir", "/path/./to/my/dir/file", "./to/my/dir/file"),
    ("/path/to/my/dir/with/file", "/path/to/other/dir/with/file", "other/dir/with/file")
  )

  private def goodPaths = Seq(

    // Normal paths, not normalized

    GoodPath(
      description = "a path",
      path = "/hello/world",
      normalize = false,
      pathAsString = "/hello/world",
      pathWithoutScheme = "/hello/world",
      parent = "/hello",
      getParent = "/hello",
      root = "/",
      name = "world",
      getFileName = "world",
      getNameCount = 2,
      isAbsolute = true,
      isDirectory = false),

    GoodPath(
      description = "a relative path",
      path = "hello/world",
      normalize = false,
      pathAsString = "hello/world",
      pathWithoutScheme = s"$pwdPathPrefix/hello/world",
      parent = pwdPathPrefix + "/hello",
      getParent = "hello",
      root = "/",
      name = "world",
      getFileName = "world",
      getNameCount = 2,
      isAbsolute = false,
      isDirectory = false),

    GoodPath(
      description = "a path with spaces",
      path = "/hello/world/with spaces",
      normalize = false,
      pathAsString = "/hello/world/with spaces",
      pathWithoutScheme = "/hello/world/with spaces",
      parent = "/hello/world",
      getParent = "/hello/world",
      root = "/",
      name = "with spaces",
      getFileName = "with spaces",
      getNameCount = 3,
      isAbsolute = true,
      isDirectory = false),

    GoodPath(
      description = "a path with encode spaces",
      path = "/hello/world/encoded%20spaces",
      normalize = false,
      pathAsString = "/hello/world/encoded%20spaces",
      pathWithoutScheme = "/hello/world/encoded%20spaces",
      parent = "/hello/world",
      getParent = "/hello/world",
      root = "/",
      name = "encoded%20spaces",
      getFileName = "encoded%20spaces",
      getNameCount = 3,
      isAbsolute = true,
      isDirectory = false),

    GoodPath(
      description = "a path with non-ascii characters",
      path = "/hello/world/with non ascii £€",
      normalize = false,
      pathAsString = "/hello/world/with non ascii £€",
      pathWithoutScheme = "/hello/world/with non ascii £€",
      parent = "/hello/world",
      getParent = "/hello/world",
      root = "/",
      name = "with non ascii £€",
      getFileName = "with non ascii £€",
      getNameCount = 3,
      isAbsolute = true,
      isDirectory = false),

    // Special paths

    GoodPath(
      description = "an empty path",
      path = "",
      normalize = false,
      pathAsString = "",
      pathWithoutScheme = pwdPath,
      parent = parentPath,
      getParent = null,
      root = "/",
      name = pwdName,
      getFileName = "",
      getNameCount = 1,
      isAbsolute = false,
      isDirectory = true),

    GoodPath(
      description = "a path from /",
      path = "/",
      normalize = false,
      pathAsString = "/",
      pathWithoutScheme = "/",
      parent = null,
      getParent = null,
      root = "/",
      name = "",
      getFileName = null,
      getNameCount = 0,
      isAbsolute = true,
      isDirectory = true),

    GoodPath(
      description = "a path from .",
      path = ".",
      normalize = false,
      pathAsString = ".",
      pathWithoutScheme = s"$pwdPathPrefix/.",
      parent = parentPath,
      getParent = null,
      root = "/",
      name = pwdName,
      getFileName = ".",
      getNameCount = 1,
      isAbsolute = false,
      isDirectory = true),

    GoodPath(
      description = "a path from ..",
      path = "..",
      normalize = false,
      pathAsString = "..",
      pathWithoutScheme = s"$pwdPathPrefix/..",
      parent = grandParentPath,
      getParent = null,
      root = "/",
      name = parentName,
      getFileName = "..",
      getNameCount = 1,
      isAbsolute = false,
      isDirectory = true),

    GoodPath(
      description = "a path including .",
      path = "/hello/world/with/./dots",
      normalize = false,
      pathAsString = "/hello/world/with/./dots",
      pathWithoutScheme = "/hello/world/with/./dots",
      parent = "/hello/world/with",
      getParent = "/hello/world/with/.",
      root = "/",
      name = "dots",
      getFileName = "dots",
      getNameCount = 5,
      isAbsolute = true,
      isDirectory = false),

    GoodPath(
      description = "a path including ..",
      path = "/hello/world/with/../dots",
      normalize = false,
      pathAsString = "/hello/world/with/../dots",
      pathWithoutScheme = "/hello/world/with/../dots",
      parent = "/hello/world",
      getParent = "/hello/world/with/..",
      root = "/",
      name = "dots",
      getFileName = "dots",
      getNameCount = 5,
      isAbsolute = true,
      isDirectory = false),

    // Normalized

    GoodPath(
      description = "an normalized empty path",
      path = "",
      normalize = true,
      pathAsString = "",
      pathWithoutScheme = pwdPath,
      parent = parentPath,
      getParent = null,
      root = "/",
      name = pwdName,
      getFileName = "",
      getNameCount = 1,
      isAbsolute = false,
      isDirectory = true),

    GoodPath(
      description = "a normalized path from /",
      path = "/",
      normalize = true,
      pathAsString = "/",
      pathWithoutScheme = "/",
      parent = null,
      getParent = null,
      root = "/",
      name = "",
      getFileName = null,
      getNameCount = 0,
      isAbsolute = true,
      isDirectory = true),

    GoodPath(
      description = "a normalized path from .",
      path = ".",
      normalize = true,
      pathAsString = "",
      pathWithoutScheme = pwdPath,
      parent = parentPath,
      getParent = null,
      root = "/",
      name = pwdName,
      getFileName = "",
      getNameCount = 1,
      isAbsolute = false,
      isDirectory = true),

    GoodPath(
      description = "a normalized path from ..",
      path = "..",
      normalize = true,
      pathAsString = "..",
      pathWithoutScheme = s"$pwdPathPrefix/..",
      parent = grandParentPath,
      getParent = null,
      root = "/",
      name = parentName,
      getFileName = "..",
      getNameCount = 1,
      isAbsolute = false,
      isDirectory = true),

    GoodPath(
      description = "a normalized path including a .",
      path = "/hello/world/with/./dots",
      normalize = true,
      pathAsString = "/hello/world/with/dots",
      pathWithoutScheme = "/hello/world/with/dots",
      parent = "/hello/world/with",
      getParent = "/hello/world/with",
      root = "/",
      name = "dots",
      getFileName = "dots",
      getNameCount = 4,
      isAbsolute = true,
      isDirectory = false),

    GoodPath(
      description = "a normalized path including ..",
      path = "/hello/world/with/../dots",
      normalize = true,
      pathAsString = "/hello/world/dots",
      pathWithoutScheme = "/hello/world/dots",
      parent = "/hello/world",
      getParent = "/hello/world",
      root = "/",
      name = "dots",
      getFileName = "dots",
      getNameCount = 3,
      isAbsolute = true,
      isDirectory = false),

    // URI

    GoodPath(
      description = "a path from a file uri",
      path = "file:///hello/world",
      normalize = false,
      pathAsString = "/hello/world",
      pathWithoutScheme = "/hello/world",
      parent = "/hello",
      getParent = "/hello",
      root = "/",
      name = "world",
      getFileName = "world",
      getNameCount = 2,
      isAbsolute = true,
      isDirectory = false),

    GoodPath(
      description = "a path from a file uri with encoded spaces",
      path = "file:///hello/world/encoded%20spaces",
      normalize = false,
      pathAsString = "/hello/world/encoded%20spaces",
      pathWithoutScheme = "/hello/world/encoded%20spaces",
      parent = "/hello/world",
      getParent = "/hello/world",
      root = "/",
      name = "encoded%20spaces",
      getFileName = "encoded%20spaces",
      getNameCount = 3,
      isAbsolute = true,
      isDirectory = false)
  )

  private def badPaths = Seq(
    BadPath("a https path", "https://hello/world", "Cannot build a local path from https://hello/world"),
    BadPath("a gcs path", "gs://hello/world", "Cannot build a local path from gs://hello/world")
  )
}
