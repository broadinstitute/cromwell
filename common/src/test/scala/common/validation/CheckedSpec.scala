package common.validation

import cats.data.NonEmptyList
import org.scalatest.{FlatSpec, Matchers}
import common.validation.Checked._

class CheckedSpec extends FlatSpec with Matchers {
  behavior of "Checked"
  
  it should "provide helper methods" in {
    5.validNelCheck shouldBe Right(5)
    "argh".invalidNelCheck[Int] shouldBe Left(NonEmptyList.one("argh"))
  }
}
