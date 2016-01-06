package cromwell.engine.io.gcs

import org.scalatest.mock.MockitoSugar
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.ExecutionContext

class NioGcsPathSpec extends FlatSpec with Matchers with MockitoSugar {

  behavior of "NioGcsPath"

  implicit val GCSFs = GcsFileSystem.getInstance(mock[GoogleCloudStorage], "gs://root")(mock[ExecutionContext])

  it should "implement toString" in {
    val absPath1 = new NioGcsPath(Array("absolute", "path", "to", "somewhere"), true)
    val relPath1 = new NioGcsPath(Array("some", "relative", "path"), false)

    absPath1.toString shouldBe "gs://absolute/path/to/somewhere"
    relPath1.toString shouldBe "some/relative/path"
  }

  it should "implement subpath" in {
    val absPath1 = new NioGcsPath(Array("absolute", "path", "to", "somewhere"), true)
    val relPath1 = new NioGcsPath(Array("some", "relative", "path"), false)

    val absSub1 = absPath1.subpath(0, 2)
    absSub1.isAbsolute shouldBe true
    absSub1.toString shouldBe "gs://absolute/path"

    val absSub2 = absPath1.subpath(1, 2)
    absSub2.isAbsolute shouldBe false
    absSub2.toString shouldBe "path"

    val relSub1 = relPath1.subpath(0, 2)
    relSub1.isAbsolute shouldBe false
    relSub1.toString shouldBe "some/relative"
  }

  it should "implement resolveSibling" in {
    val absPath1 = new NioGcsPath(Array("absolute", "path", "to", "somewhere"), true)
    val absPath2 = new NioGcsPath(Array("absolute", "location"), true)
    val relPath1 = new NioGcsPath(Array("some", "relative", "path"), false)
    val relPath2 = new NioGcsPath(Array("another", "relative", "resource", "path"), false)

    val absSibling = absPath1.resolveSibling("somewhere else")
    absSibling.isAbsolute shouldBe true
    absSibling.toString shouldBe "gs://absolute/path/to/somewhere else"

    val absSiblingPath = absPath1.resolveSibling(relPath1)
    absSiblingPath.isAbsolute shouldBe true
    absSiblingPath.toString shouldBe "gs://absolute/path/to/some/relative/path"

    val absRel = relPath1.resolveSibling("other path")
    absRel.isAbsolute shouldBe false
    absRel.toString shouldBe "some/relative/other path"

    val absRelPath = relPath1.resolveSibling(relPath2)
    absRelPath.isAbsolute shouldBe false
    absRelPath.toString shouldBe "some/relative/another/relative/resource/path"
  }

  it should "implement resolve" in {
    val absPath1 = new NioGcsPath(Array("absolute", "path", "to", "somewhere"), true)
    val absPath2 = new NioGcsPath(Array("absolute", "location"), true)
    val relPath1 = new NioGcsPath(Array("some", "relative", "path"), false)
    val relPath2 = new NioGcsPath(Array("another", "relative", "resource", "path"), false)

    val absToRel = absPath1.resolve(relPath1)
    absToRel.isAbsolute shouldBe true
    absToRel.toString shouldBe "gs://absolute/path/to/somewhere/some/relative/path"

    val absToAbs = absPath1.resolve(absPath2)
    absToAbs.isAbsolute shouldBe true
    absToAbs.toString shouldBe "gs://absolute/location"

    val relToAbs = relPath1.resolve(absPath1)
    relToAbs.isAbsolute shouldBe true
    relToAbs.toString shouldBe "gs://absolute/path/to/somewhere"

    val relToRel = relPath1.resolve(relPath2)
    relToRel.isAbsolute shouldBe false
    relToRel.toString shouldBe "some/relative/path/another/relative/resource/path"
  }

}
