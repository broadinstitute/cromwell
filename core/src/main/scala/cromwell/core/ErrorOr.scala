package cromwell.core

import cats.data.Validated.{Invalid, Valid}
import cats.data.{NonEmptyList, Validated}

object ErrorOr {
  type ErrorOr[A] = Validated[NonEmptyList[String], A]

  implicit class ShortCircuitingFlatMap[A](val fa: ErrorOr[A]) extends AnyVal {
    /**
      * Not consistent with `Applicative#ap` but useful in for comprehensions.
      *
      * @see http://typelevel.org/cats/tut/validated.html#of-flatmaps-and-xors
      */
    def flatMap[B](f: A => ErrorOr[B]): ErrorOr[B] = {
      fa match {
        case Valid(v) => f(v)
        case i @ Invalid(_) => i
      }
    }
  }
}
