package cromwell.filesystems.oss.nio

import cromwell.core.TestKitSuite
import org.scalatest.TryValues._

import scala.util.{Failure, Success, Try}

case class ValidPath(pathAsString: String,
                     parent: Option[String],
                     fileName: Option[String],
                     nameCount: Int,
                     root: Option[String] = Some(UnixPath.ROOT_PATH.toString),
                     isRoot: Boolean = false,
                     isAbsolute: Boolean = true,
                     hasTrailingSeparator: Boolean = false,
                     likeDir: Boolean = false,
                    )

case class SubPath(pathAsString: String,
                   beginIndex: Int,
                   endIndex: Int,
                   subPath: Try[String],
                   description: String = ""
                  )

case class ResolvePath(pathAsString: String,
                       other: String,
                       resolved: String,
                       description: String = ""
                      )

case class ResolveSiblingPath(pathAsString: String,
                              other: String,
                              resolved: String,
                              description: String = ""
                             )

case class NormalizePath(pathAsString: String,
                         normalized: String
                        )

case class RelativizePath(pathAsString: String,
                          other: String,
                          relative: String
                         )

class UnixPathSpec extends TestKitSuite with OssNioUtilSpec  {

  validPaths foreach { path =>
    it should behave like verifyValidPath(path)
  }

  subPaths foreach { path =>
    it should behave like verifySubPath(path)
  }

  resolvePaths foreach { path =>
    it should behave like verifyResolvePath(path)
  }

  resolveSiblingPaths foreach { path =>
    it should behave like verifyResolveSiblingPath(path)
  }

  normalizePaths foreach { path =>
    it should behave like verifyNormalizePath(path)
  }

  relativizePaths foreach { path =>
    it should behave like verifyRelativePath(path)
  }


  def verifyValidPath(path: ValidPath) = {
    behavior of s"Verify a valid UnixPath ${path.pathAsString}"

    val clue = s"pathAsString: ${path.pathAsString}"

    val unixPath = UnixPath(path.pathAsString)

    it should "match expected parent" in
      withClue(clue) {
        unixPath.getParent map {_.toString} shouldEqual path.parent
      }

    it should "match expected name count" in
      withClue(clue) {
        unixPath.getNameCount shouldEqual path.nameCount
      }

    it should "match expected file name" in
      withClue(clue) {
        unixPath.getFileName map {_.toString} shouldEqual path.fileName
      }

    it should "match expected root" in
      withClue(clue) {
        unixPath.getRoot map {_.toString} shouldEqual path.root
      }

    it should "match expected isRoot" in
      withClue(clue) {
        unixPath.isRoot shouldBe path.isRoot
      }

    it should "match expected isAbsolute" in
      withClue(clue) {
        unixPath.isAbsolute shouldBe path.isAbsolute
      }

    it should "match expected hasTrailingSeparator" in
      withClue(clue) {
        unixPath.hasTrailingSeparator shouldBe path.hasTrailingSeparator
      }

    it should "match expected seemsLikeDirectory" in
      withClue(clue) {
        unixPath.seemsLikeDirectory() shouldBe path.likeDir
      }
  }

  def verifySubPath(path: SubPath) = {
    val clue = s"path ${path.pathAsString} beginIndex ${path.beginIndex} endIndex ${path.endIndex}"

    behavior of s"Verify a unix path's sub path ${clue}"

    val unixPath  = UnixPath(path.pathAsString)

    it should "match sub path" in
      withClue(clue) {
        val maybeRes = unixPath.subPath(path.beginIndex, path.endIndex)
        path.subPath match {
          case Success(sub: String) =>
            maybeRes.success.value.toString shouldEqual(sub)
          case Failure(_) =>
            maybeRes.failure
        }
      }
  }

  def verifyResolvePath(path: ResolvePath) = {
    val clue = s"path ${path.pathAsString} other ${path.other}"

    behavior of s"Verify resolving a path on a UnixPath ${clue}"

    val me  = UnixPath(path.pathAsString)
    val other = UnixPath(path.other)

    it should "match resolved path" in
      withClue(clue) {
        me.resolve(other).toString shouldEqual(path.resolved)
      }
  }

  def verifyResolveSiblingPath(path: ResolveSiblingPath) = {
    val clue = s"path ${path.pathAsString} other ${path.other}"

    behavior of s"Verify resolving sibling on a UnixPath ${clue}"

    val me  = UnixPath(path.pathAsString)
    val other = UnixPath(path.other)

    it should "match expected sibling path" in
      withClue(clue) {
        me.resolveSibling(other).toString shouldEqual(path.resolved)
      }
  }


  def verifyNormalizePath(path: NormalizePath) = {
    val clue = s"path ${path.pathAsString}"

    behavior of s"Verify normalize a UnixPath ${clue}"

    val me  = UnixPath(path.pathAsString)

    it should "match expected normalized path" in
      withClue(clue) {
        me.normalize().toString shouldEqual(path.normalized)
      }
  }


  def verifyRelativePath(path: RelativizePath) = {
    val clue = s"path ${path.pathAsString} other ${path.other}"

    behavior of s"Verify resolving relativize on a UnixPath ${clue}"

    val me  = UnixPath(path.pathAsString)
    val other = UnixPath(path.other)

    it should "match resolved path" in
      withClue(clue) {
        me.relativize(other).toString shouldEqual(path.relative)
      }
  }

  private[this] def validPaths = Seq(
    ValidPath(
      pathAsString = "/bcs-dir/bcs-file",
      parent = Some("/bcs-dir/"),
      fileName = Some("bcs-file"),
      nameCount = 2
    ),
    ValidPath(
      pathAsString = "/bcs-dir/bcs-dir1/",
      parent = Some("/bcs-dir/"),
      fileName = Some("bcs-dir1"),
      nameCount = 2,
      hasTrailingSeparator = true,
      likeDir = true
    ),
    ValidPath(
      pathAsString = "bcs-file",
      parent = None,
      fileName = Some("bcs-file"),
      nameCount = 1,
      root = None,
      isAbsolute = false,
    )
  )

  private[this] def subPaths = Seq(
    SubPath(
      pathAsString = "/bcs-dir/bcs-dir1/bcs-dir2",
      beginIndex = 0,
      endIndex = 1,
      subPath = Success("bcs-dir"),
      description = "valid slice"
    ),
    SubPath(
      pathAsString = "/bcs-dir/bcs-dir1/bcs-dir2",
      beginIndex = 1,
      endIndex = 0,
      subPath = Failure(new IllegalArgumentException()),
      description = "invalid index"
    ),
    SubPath(
      pathAsString = "/bcs-dir/bcs-dir1/bcs-dir2",
      beginIndex = 1,
      endIndex = 10,
      subPath = Success("bcs-dir1/bcs-dir2"),
      description = "valid index"
    ),
    SubPath(
      pathAsString = "/bcs-dir/bcs-dir1/bcs-dir2",
      beginIndex = 3,
      endIndex = 10,
      subPath = Success(""),
      description = "valid index"
    )
  )

  private[this] def resolvePaths = Seq(
    ResolvePath(
      pathAsString = "/bcs-dir/bcs-dir1",
      other = "",
      resolved = "/bcs-dir/bcs-dir1"
    ),
    ResolvePath(
      pathAsString = "/bcs-dir/bcs-dir1",
      other = "/bcs-dir2/bcs-dir3",
      resolved = "/bcs-dir2/bcs-dir3"
    ),
    ResolvePath(
      pathAsString = "/bcs-dir/bcs-dir1/",
      other = "bcs-file",
      resolved = "/bcs-dir/bcs-dir1/bcs-file"
    ),
    ResolvePath(
      pathAsString = "/bcs-dir/bcs-dir1",
      other = "bcs-file",
      resolved = "/bcs-dir/bcs-dir1/bcs-file"
    ),
  )

  private[this] def resolveSiblingPaths = Seq(
    ResolveSiblingPath(
      pathAsString = "/bcs-dir/bcs-file1",
      other = "bcs-file2",
      resolved = "/bcs-dir/bcs-file2"
    ),
    ResolveSiblingPath(
      pathAsString = "/",
      other = "bcs-file2",
      resolved = "bcs-file2"
    ),
    ResolveSiblingPath(
      pathAsString = "",
      other = "bcs-file2",
      resolved = "bcs-file2"
    ),
    ResolveSiblingPath(
      pathAsString = "/bcs-dir/bcs-file1",
      other = "/bcs-file2",
      resolved = "/bcs-file2"
    ),
  )

  private[this] def normalizePaths = Seq(
    NormalizePath(
      "/bcs-dir/.",
      "/bcs-dir/"
    ),
    NormalizePath(
      "/bcs-dir/../bcs-file",
      "/bcs-file"
    ),
    NormalizePath(
      "/bcs-dir/./bcs-file",
      "/bcs-dir/bcs-file"
    ),
    NormalizePath(
      "/bcs-dir/./bcs-dir1/",
      "/bcs-dir/bcs-dir1/"
    ),
    NormalizePath(
      "../bcs-dir/bcs-dir1/",
      "bcs-dir/bcs-dir1/"
    )
  )

  private[this] def relativizePaths = Seq(
    RelativizePath(
      "/bcs-dir1/bcs-file1",
      "/bcs-dir1/bcs-file2",
      "../bcs-file2"
    ),
    RelativizePath(
      "/bcs-dir1/bcs-file1",
      "/bcs-file2",
      "../../bcs-file2"
    )
  )
}
