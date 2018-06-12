package cromwell.core.path

import org.scalatest.Matchers._
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop._
import org.scalatest.{FlatSpecLike, Tag}

case class GoodPath(description: String,
                    path: String,
                    normalize: Boolean,
                    pathAsString: String,
                    pathWithoutScheme: String,
                    parent: String,
                    getParent: String,
                    root: String,
                    name: String,
                    getFileName: String,
                    getNameCount: Int,
                    isAbsolute: Boolean,
                    isDirectory: Boolean)

case class BadPath(description: String, path: String, exceptionMessage: String)

/**
  * Tests various methods of Path objects.
  */
trait PathBuilderSpecUtils {
  this: FlatSpecLike =>

  def truncateCommonRoots(builder: => PathBuilder,
                          pathsToTruncate: TableFor3[String, String, String],
                          tag: Tag = PathBuilderSpecUtils.PathTest): Unit = {
    behavior of s"PathCopier"

    it should "truncate common roots" taggedAs tag in {
      forAll(pathsToTruncate) { (context, file, relative) =>
        val contextPath = builder.build(context).get
        val filePath = builder.build(file).get
        val actual = PathCopier.truncateCommonRoot(contextPath, filePath)
        actual should be(relative)
      }
    }
  }

  def buildGoodPath(builder: => PathBuilder, goodPath: GoodPath, tag: Tag = PathBuilderSpecUtils.PathTest): Unit = {
    behavior of s"Building ${goodPath.description}"

    lazy val path = {
      val path = builder.build(goodPath.path).get
      if (goodPath.normalize) path.normalize() else path
    }

    it should "match expected pathAsString" taggedAs tag in
      withClue(s"for path ${goodPath.path}${if (goodPath.normalize) " (normalized)" else ""}:") {
        path.pathAsString should be(goodPath.pathAsString)
      }

    it should "match expected pathWithoutScheme" taggedAs tag in
      withClue(s"for path ${goodPath.path}${if (goodPath.normalize) " (normalized)" else ""}:") {
        path.pathWithoutScheme should be(goodPath.pathWithoutScheme)
      }

    it should "match expected parent" taggedAs tag in
      withClue(s"for path ${goodPath.path}${if (goodPath.normalize) " (normalized)" else ""}:") {
        Option(path.parent).map(_.pathAsString).orNull should be(goodPath.parent)
      }

    it should "match expected getParent" taggedAs tag in
      withClue(s"for path ${goodPath.path}${if (goodPath.normalize) " (normalized)" else ""}:") {
        Option(path.getParent).map(_.pathAsString).orNull should be(goodPath.getParent)
      }

    it should "match expected root" taggedAs tag in
      withClue(s"for path ${goodPath.path}${if (goodPath.normalize) " (normalized)" else ""}:") {
        Option(path.root).map(_.pathAsString).orNull should be(goodPath.root)
      }

    it should "match expected name" taggedAs tag in
      withClue(s"for path ${goodPath.path}${if (goodPath.normalize) " (normalized)" else ""}:") {
        path.name should be(goodPath.name)
      }

    it should "match expected getFileName" taggedAs tag in
      withClue(s"for path ${goodPath.path}${if (goodPath.normalize) " (normalized)" else ""}:") {
        Option(path.getFileName).map(_.pathAsString).orNull should be(goodPath.getFileName)
      }

    it should "match expected getNameCount" taggedAs tag in
      withClue(s"for path ${goodPath.path}${if (goodPath.normalize) " (normalized)" else ""}:") {
        path.getNameCount should be(goodPath.getNameCount)
      }

    it should "match expected isAbsolute" taggedAs tag in
      withClue(s"for path ${goodPath.path}${if (goodPath.normalize) " (normalized)" else ""}:") {
        path.isAbsolute should be(goodPath.isAbsolute)
      }

    it should "match expected isDirectory" taggedAs tag in
      withClue(s"for path ${goodPath.path}${if (goodPath.normalize) " (normalized)" else ""}:") {
        path.isDirectory should be(goodPath.isDirectory)
      }
  }

  def buildBadPath(builder: => PathBuilder, badPath: BadPath, tag: Tag = PathBuilderSpecUtils.PathTest): Unit = {
    behavior of s"Building ${badPath.description}"

    it should "fail to build a path" taggedAs tag in
      withClue(s"for path ${badPath.path}:") {
        val exception = intercept[Exception](builder.build(badPath.path).get)
        exception.getMessage should be(badPath.exceptionMessage)
      }
  }

}

object PathBuilderSpecUtils {
  val PathTest = Tag("PathTest")
  val AwsTest = Tag("AwsTest")
}
