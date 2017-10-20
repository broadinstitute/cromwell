package common.validation

import cats.data.NonEmptyList
import cats.syntax.either._
import common.Checked

object Checked {
  implicit class ValidCheck[A](val obj: A) extends AnyVal {
    def validNelCheck: Checked[A] = obj.asRight[NonEmptyList[String]]
  }

  implicit class InvalidCheck(val obj: String) extends AnyVal {
    def invalidNelCheck[A]: Checked[A] = NonEmptyList.one(obj).asLeft[A]
  }
}
