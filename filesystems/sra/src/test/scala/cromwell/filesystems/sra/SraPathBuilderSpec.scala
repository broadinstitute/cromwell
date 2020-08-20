package cromwell.filesystems.sra

import cromwell.core.TestKitSuite
import cromwell.core.path._
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

class SraPathBuilderSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with PathBuilderSpecUtils {
  behavior of "SraPathBuilder"

  goodPaths foreach {
    good => it should behave like buildGood(good)
  }

  badPaths foreach {
    bad => it should behave like buildBad(bad)
  }

  private def buildGood(good: Good): Unit = {
    behavior of s"Building ${good.description}"

    val built = pathBuilder.build(good.in)
    it should "successfully construct an SraPath" in {
      assert(built.isSuccess)
    }

    val path = built.get
    val sraPath = path.asInstanceOf[SraPath]

    it should "match expected accession" in {
      sraPath.accession should be (good.accession)
    }

    it should "match expected path" in {
      sraPath.path should be (good.path)
    }
  }

  private def buildBad(bad: Bad): Unit = {
    behavior of s"Building ${bad.description}"

    it should "fail to build an SraPath" in {
      pathBuilder.build(bad.in).isSuccess should be (false)
    }
  }

  private lazy val pathBuilder = new SraPathBuilder

  private case class Good(description: String,
                          in: String,
                          accession: String,
                          path: String)
  private def goodPaths = Seq(
    Good(
      description = "well-formed path",
      in = "sra://SXP000001/asdf.bam",
      accession = "SXP000001",
      path = "asdf.bam",
    ),
    Good(
      description = "nested path",
      in = "sra://SRA42424242/first/second/third.bam.bai",
      accession = "SRA42424242",
      path = "first/second/third.bam.bai",
    ),
    Good(
      description = "path with spaces",
      in = "sra://SXP111111/top level/nested level/file name.bz2",
      accession = "SXP111111",
      path = "top level/nested level/file name.bz2",
    ),
  )

  private case class Bad(description: String, in: String)
  private def badPaths = Seq(
    Bad(
      description = "not an SRA path",
      in = "gcs://some/gcs/path/thing.txt",
    ),
    Bad(
      description = "missing accession",
      in = "sra://",
    ),
    Bad(
      description = "missing path within accession",
      in = "sra://SRA00001",
    ),
  )
}
