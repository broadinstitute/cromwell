package common.assertion

import cats.data.Validated.{Invalid, Valid}
import common.Checked
import common.validation.ErrorOr.ErrorOr
import org.scalatest.{Assertion, Matchers}

object ErrorOrAssertions {
  implicit class ErrorOrWithAssertions[A](errorOr: ErrorOr[A]) extends Matchers {
    def shouldBeValid(other: A): Assertion = errorOr match {
      case Valid(mine) => mine shouldBe other
      case Invalid(e) => fail(s"Invalid ErrorOr: ${e.toList.mkString(", ")}")
    }

    def shouldBeInvalid(errors: String*): Assertion = errorOr match {
      case Valid(mine) =>
        fail(s"Unexpectedly 'Valid' Checked: $mine")
      case Invalid(e) =>
        e.toList should be(errors.toList)
    }
  }

  implicit class CheckedWithAssertions[A](checked: Checked[A]) extends Matchers {
    def shouldBeValid(other: A): Assertion = checked match {
      case Right(mine) => mine shouldBe other
      case Left(e) => fail(s"Invalid ErrorOr: ${e.toList.mkString(", ")}")
    }

    def shouldBeValidPF(f: PartialFunction[A, Assertion]): Assertion = checked match {
      case Right(mine) if f.isDefinedAt(mine) => f(mine)
      case Right(other) => fail(s"Unexpected value in a valid Checked: $other")
      case Left(e) =>
        fail(s"Invalid ErrorOr: ${e.toList.mkString(", ")}")
    }

    def shouldBeInvalid(errors: String*): Assertion = checked match {
      case Right(mine) =>
        fail(s"Unexpectedly 'Valid' Checked: $mine")
      case Left(e) =>
        e.toList should be(errors.toList)
    }
  }
}
