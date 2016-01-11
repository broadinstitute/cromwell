package cromwell.engine.io.gcs

import java.nio.file.Path

import org.scalatest.mock.MockitoSugar
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.ExecutionContext

class NioGcsPathSpec extends FlatSpec with Matchers with MockitoSugar {

  behavior of "NioGcsPath"

  implicit val GCSFs = GcsFileSystem.instance(mock[GoogleCloudStorage], "gs://absolute")(mock[ExecutionContext])

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

  it should "implement getName" in {
    val absPath1 = new NioGcsPath(Array("absolute", "path", "to", "somewhere"), true)
    val relPath1 = new NioGcsPath(Array("some", "relative", "path"), false)

    val nameAbs1 = absPath1.getName(0)
    nameAbs1.isAbsolute shouldBe true
    nameAbs1.toString shouldBe "gs://absolute"

    val nameAbs2 = absPath1.getName(1)
    nameAbs2.isAbsolute shouldBe false
    nameAbs2.toString shouldBe "path"

    val nameRel1 = relPath1.getName(0)
    nameRel1.isAbsolute shouldBe false
    nameRel1.toString shouldBe "some"

    val nameRel2 = relPath1.getName(1)
    nameRel2.isAbsolute shouldBe false
    nameRel2.toString shouldBe "relative"
  }

  it should "implement getParent" in {
    val empty = new NioGcsPath(Array.empty[String], true)
    val singleton = new NioGcsPath(Array("singleton"), true)
    val absPath1 = new NioGcsPath(Array("absolute", "path", "to", "somewhere"), true)
    val relPath1 = new NioGcsPath(Array("some", "relative", "path"), false)

    val parentAbs1 = absPath1.getParent
    parentAbs1.isAbsolute shouldBe true
    parentAbs1.toString shouldBe "gs://absolute/path/to"

    empty.getParent shouldBe null
    singleton.getParent shouldBe null

    val nameRel1 = relPath1.getParent
    nameRel1.isAbsolute shouldBe false
    nameRel1.toString shouldBe "some/relative"
  }

  it should "implement toAbsolutePath" in {
    val absPath1 = new NioGcsPath(Array("absolute", "path", "to", "somewhere"), true)
    val relPath1 = new NioGcsPath(Array("some", "relative", "path"), false)

    val abs = absPath1.toAbsolutePath
    abs.isAbsolute shouldBe true
    abs.toString shouldBe "gs://absolute/path/to/somewhere"

    val rel = relPath1.toAbsolutePath
    rel.isAbsolute shouldBe true
    rel.toString shouldBe "gs://absolute/some/relative/path"
  }

  it should "implement getNameCount" in {
    val empty = new NioGcsPath(Array.empty[String], true)
    val singleton = new NioGcsPath(Array("singleton"), true)
    val absPath1 = new NioGcsPath(Array("absolute", "path", "to", "somewhere"), true)
    val relPath1 = new NioGcsPath(Array("some", "relative", "path"), false)

    absPath1.getNameCount shouldBe 4
    relPath1.getNameCount shouldBe 3
    empty.getNameCount shouldBe 0
    singleton.getNameCount shouldBe 1
  }

  it should "implement getFileName" in {
    val empty = new NioGcsPath(Array.empty[String], true)
    val singletonAbs = new NioGcsPath(Array("singleton"), true)
    val singletonRel = new NioGcsPath(Array("singleton"), false)
    val absPath1 = new NioGcsPath(Array("absolute", "path", "to", "somewhere"), true)
    val relPath1 = new NioGcsPath(Array("some", "relative", "path"), false)

    val emptyFileName = empty.getFileName
    emptyFileName shouldBe null

    val singletonAbsFileName = singletonAbs.getFileName
    singletonAbsFileName.isAbsolute shouldBe true
    singletonAbsFileName.toString shouldBe "gs://singleton"

    val singletonRelFileName = singletonRel.getFileName
    singletonRelFileName.isAbsolute shouldBe false
    singletonRelFileName.toString shouldBe "singleton"

    val relFileName = relPath1.getFileName
    relFileName.isAbsolute shouldBe false
    relFileName.toString shouldBe "path"

    val absFileName = absPath1.getFileName
    absFileName.isAbsolute shouldBe false
    absFileName.toString shouldBe "somewhere"
  }

  it should "implement getRoot" in {
    val empty = new NioGcsPath(Array.empty[String], true)
    val singletonAbs = new NioGcsPath(Array("singleton"), true)
    val singletonRel = new NioGcsPath(Array("singleton"), false)
    val absPath1 = new NioGcsPath(Array("absolute", "path", "to", "somewhere"), true)
    val relPath1 = new NioGcsPath(Array("some", "relative", "path"), false)

    an[IllegalStateException] shouldBe thrownBy(empty.getRoot)

    val singletonAbsFileName = singletonAbs.getRoot
    singletonAbsFileName.isAbsolute shouldBe true
    singletonAbsFileName.toString shouldBe "gs://singleton"

    val singletonRelFileName = singletonRel.getRoot
    singletonRelFileName.isAbsolute shouldBe true
    singletonRelFileName.toString shouldBe "gs://absolute"

    val relFileName = relPath1.getRoot
    relFileName.isAbsolute shouldBe true
    relFileName.toString shouldBe "gs://absolute"

    val absFileName = absPath1.getRoot
    absFileName.isAbsolute shouldBe true
    absFileName.toString shouldBe "gs://absolute"
  }

  it should "implement getIterator" in {
    val empty = new NioGcsPath(Array.empty[String], true)
    val singletonAbs = new NioGcsPath(Array("singleton"), true)
    val singletonRel = new NioGcsPath(Array("singleton"), false)
    val absPath1 = new NioGcsPath(Array("absolute", "path", "to", "somewhere"), true)
    val relPath1 = new NioGcsPath(Array("some", "relative", "path"), false)

    empty.iterator().hasNext shouldBe false

    val singletonAbsIterator = singletonAbs.iterator()
    val nextAbsSingleton: Path = singletonAbsIterator.next()
    nextAbsSingleton.isAbsolute shouldBe false
    nextAbsSingleton.toString shouldBe "singleton"
    singletonAbsIterator.hasNext shouldBe false

    val singletonRelIterator = singletonRel.iterator()
    val nextRelSingleton: Path = singletonRelIterator.next()
    nextRelSingleton.isAbsolute shouldBe false
    nextRelSingleton.toString shouldBe "singleton"
    singletonRelIterator.hasNext shouldBe false

    val relIterator = relPath1.iterator()
    val nextRel: Path = relIterator.next()
    nextRel.isAbsolute shouldBe false
    nextRel.toString shouldBe "some"
    relIterator.next().toString shouldBe "relative"
    relIterator.next().toString shouldBe "path"
    relIterator.hasNext shouldBe false

    val absIterator = absPath1.iterator()
    val absRel: Path = absIterator.next()
    absRel.isAbsolute shouldBe false
    absRel.toString shouldBe "absolute"
    absIterator.next().toString shouldBe "path"
    absIterator.next().toString shouldBe "to"
    absIterator.next().toString shouldBe "somewhere"
    absIterator.hasNext shouldBe false
  }

  it should "implement startsWith" in {
    val empty = new NioGcsPath(Array.empty[String], false)
    val singletonAbs = new NioGcsPath(Array("absolute"), true)
    val singletonRel = new NioGcsPath(Array("absolute"), false)

    val absPath = new NioGcsPath(Array("absolute", "path", "to", "somewhere"), true)
    val startsWithAbsPath = new NioGcsPath(Array("absolute", "path", "to"), true)
    val doesntStartsWithAbsPath = new NioGcsPath(Array("absolute", "path", "to", "another", "place"), true)
    val absPathStartingLikeRel = new NioGcsPath(Array("some", "relative", "path"), true)

    val relPath = new NioGcsPath(Array("some", "relative", "path"), false)
    val startsWithRelPath = new NioGcsPath(Array("some", "relative"), false)
    val doesntStartsWithRelPath = new NioGcsPath(Array("some", "relative", "other", "path"), false)
    val relPathStartingLikeAbs = new NioGcsPath(Array("absolute", "path", "to"), false)

    val paths = Table(
      ("path1", "path2", "result"),
      (empty, empty, true),
      (empty, absPath, false),
      (singletonAbs, singletonAbs, true),
      (absPath, startsWithAbsPath, true),
      (absPath, doesntStartsWithAbsPath, false),
      (absPath, relPathStartingLikeAbs, true),
      (absPath, relPath, false),
      (relPath, startsWithRelPath, true),
      (relPath, doesntStartsWithRelPath, false),
      (relPath, absPathStartingLikeRel, false),
      (relPath, absPath, false)
    )

    forAll(paths) { (p1, p2, res) =>
      val startsWith: Boolean = p1.startsWith(p2)
      startsWith shouldBe res
      val startsWith1: Boolean = p1.startsWith(p2.toString)
      startsWith1 shouldBe res
    }
  }

  it should "implement endsWith" in {
    val empty = new NioGcsPath(Array.empty[String], false)
    val singletonAbs = new NioGcsPath(Array("absolute"), true)
    val singletonRel = new NioGcsPath(Array("absolute"), false)

    val absPath = new NioGcsPath(Array("absolute", "path", "to", "somewhere"), true)
    val doesntEndWithAbsPath = new NioGcsPath(Array("absolute", "path", "to", "another", "place"), true)
    val absPathEndingLikeRel = new NioGcsPath(Array("relative", "path"), true)

    val relPath = new NioGcsPath(Array("some", "relative", "path"), false)
    val endsWithRelPath = new NioGcsPath(Array("relative", "path"), false)
    val doesntStartsWithRelPath = new NioGcsPath(Array("relative", "other", "path"), false)
    val relPathEndingLikeAbs = new NioGcsPath(Array("path", "to", "somewhere"), false)

    val paths = Table(
      ("path1", "path2", "result"),
      (empty, empty, true),
      (empty, absPath, false),
      (singletonAbs, singletonAbs, true),
      (absPath, absPath, true),
      (absPath, doesntEndWithAbsPath, false),
      (absPath, relPathEndingLikeAbs, true),
      (absPath, relPath, false),
      (relPath, endsWithRelPath, true),
      (relPath, doesntStartsWithRelPath, false),
      (relPath, absPathEndingLikeRel, false),
      (relPath, absPath, false)
    )

    forAll(paths) { (p1, p2, res) =>
      p1.endsWith(p2) shouldBe res
      p1.endsWith(p2.toString) shouldBe res
    }
  }

}
