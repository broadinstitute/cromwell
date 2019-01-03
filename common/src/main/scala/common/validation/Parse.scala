package common.validation

import cats.data.EitherT.fromEither
import cats.data.{EitherT, NonEmptyList, ValidatedNel}
import cats.effect.IO
import common.Checked
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._

import scala.util.{Failure, Success, Try}

object Parse {

  type Parse[A] = EitherT[IO, NonEmptyList[String], A]

  type ParseValidated[A] = IO[ValidatedNel[String, A]]

  def error[A](error: String, tail: String*): Parse[A] = EitherT.leftT {
    NonEmptyList.of(error, tail: _*)
  }
  
  def pure[A](value: A): Parse[A] = EitherT.pure(value)

  // If anyone has a magic import that does exactly what these helpers do, please replace, thx!

  def errorOrParse[A](f: => ErrorOr[A]): Parse[A] = fromEither[IO](f.toEither)

  def checkedParse[A](f: => Checked[A]): Parse[A] = fromEither[IO](f)

  def tryParse[A](f: => Try[A]): Parse[A] = errorOrParse(f.toErrorOr)

  def goParse[A](f: => A): Parse[A] = tryParse(Try(f))

  implicit class ValidParse[A](val obj: A) extends AnyVal {
    def validParse: Parse[A] = EitherT.pure(obj)
  }

  implicit class InvalidParse(val obj: String) extends AnyVal {
    def invalidParse[A]: Parse[A] = error(obj)
  }

  implicit class EnhancedParse[A](val p: Parse[A]) extends AnyVal {
    import cats.syntax.either._
    def toChecked: Checked[A] = {
      Try(p.value.unsafeRunSync()) match {
        case Success(r) => r
        case Failure(f) => NonEmptyList.one(f.getMessage).asLeft
      }
    }
    def toErrorOr: ErrorOr[A] = toChecked.toValidated
  }
}
