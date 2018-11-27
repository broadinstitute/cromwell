package cwl

import cats.data.NonEmptyList
import cwl.CwlDecoder._
import org.scalatest.{FlatSpec, Matchers}

class CwlDecoderSpec extends FlatSpec with Matchers {

  behavior of "cwl decoder"

  import TestSetup._

  it should "read nested workflow" in {
      decodeCwlFile(rootPath / "nestedworkflows.cwl").
        value.
        unsafeRunSync match {
        case Right(cwl) =>
          val wf = cwl.select[Workflow].get
          wf.steps.flatMap(_.run.select[String].toList).length shouldBe 0
        case Left(other) => fail(other.toList.mkString(", "))
      }
  }

  it should "fail to parse broken linked cwl" in {
    val result: Either[NonEmptyList[String], Cwl] = decodeCwlFile(rootPath / "brokenlinks.cwl").
      value.
      unsafeRunSync

    result.isLeft shouldBe true
    withClue(s"Actual message: ${result.left.get.head}") {
      "Field `run` contains undefined reference to `file://.+/wrong.cwl`".r.findAllIn(result.left.get.head).size shouldBe 1
      "Field `run` contains undefined reference to `file://.+/wrong2.cwl`".r.findAllIn(result.left.get.head).size shouldBe 1
    }
  }

  it should "fail to parse invalid linked cwl" in {
      decodeCwlFile(rootPath/"links_dont_parse.cwl").
        value.
        unsafeRunSync match {
        case Left(errors) =>
          errors.filter(_.contains("bad.cwl")).size + errors.filter(_.contains("bad2.cwl")).size shouldBe 1
        case Right(_) => fail("should not have passed!")
      }
  }
}
