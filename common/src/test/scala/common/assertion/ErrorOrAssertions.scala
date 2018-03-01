package common.assertion

import cats.data.Validated.{Invalid, Valid}
import common.validation.ErrorOr.ErrorOr
import org.scalatest.{Assertion, Matchers}

object ErrorOrAssertions {
  implicit class ErrorOrWithAssertions[A](e: ErrorOr[A]) extends Matchers {
    def shouldBeValid(other: A): Assertion = e match {
      case Valid(mine) => mine shouldBe other
      case Invalid(e) => fail(s"Invalid ErrorOr: ${e.toList.mkString(", ")}")
    }
  }
}
