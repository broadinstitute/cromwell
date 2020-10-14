package common.validation

import cats.data.NonEmptyList
import common.assertion.CromwellTimeoutSpec
import common.validation.Checked._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class CheckedSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  behavior of "Checked"
  
  it should "provide helper methods" in {
    5.validNelCheck shouldBe Right(5)
    "argh".invalidNelCheck[Int] shouldBe Left(NonEmptyList.one("argh"))
  }
}
