package cwl

import cwl.CwlDecoder._
import cwl.TestSetup._
import org.scalatest._


class ParseBigThreeSpec extends FlatSpec with Matchers {
  val namespace = "cwl"

  it should "parse 1st tool" in {

  decodeAllCwl(rootPath/"1st-tool.cwl").
    value.
    unsafeRunSync.
    isRight shouldBe true
  }

  it should "parse first workflow" in {
    decodeAllCwl(rootPath/"1st-workflow.cwl").
      value.
      unsafeRunSync.
      isRight shouldBe true
  }

  it should "parse env cwl" in {
    decodeAllCwl(rootPath/"env.cwl").
      value.
      unsafeRunSync.
      isRight shouldBe true
  }
}
