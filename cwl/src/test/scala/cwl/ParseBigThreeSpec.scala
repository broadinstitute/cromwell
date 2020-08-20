package cwl

import cwl.CwlDecoder._
import cwl.TestSetup._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class ParseBigThreeSpec extends AnyFlatSpec with Matchers {
  val namespace = "cwl"

  it should "parse 1st tool" in {

  decodeCwlFile(rootPath/"1st-tool.cwl").
    value.
    unsafeRunSync.
    isRight shouldBe true
  }

  it should "parse first workflow" in {
    decodeCwlFile(rootPath/"1st-workflow.cwl").
      value.
      unsafeRunSync.
      isRight shouldBe true
  }

  it should "parse env cwl" in {
    decodeCwlFile(rootPath/"env.cwl").
      value.
      unsafeRunSync.
      isRight shouldBe true
  }
}
