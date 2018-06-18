package common.assertion

import cats.data.Validated.{Invalid, Valid}
import common.Checked
import common.validation.ErrorOr.ErrorOr
import org.scalatest.{Assertion, Matchers}

object ErrorOrAssertions {
  implicit class ErrorOrWithAssertions[A](e: ErrorOr[A]) extends Matchers {
    def shouldBeValid(other: A): Assertion = e match {
      case Valid(mine) => mine shouldBe other
      case Invalid(e) => fail(s"Invalid ErrorOr: ${e.toList.mkString(", ")}")
    }
  }

  implicit class CheckedWithAssertions[A](e: Checked[A]) extends Matchers {
    def shouldBeValid(other: A): Assertion = e match {
      case Right(mine) => mine shouldBe other
      case Left(e) => fail(s"Invalid ErrorOr: ${e.toList.mkString(", ")}")
    }
  }
}
