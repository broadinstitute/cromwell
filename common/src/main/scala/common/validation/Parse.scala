package common.validation

import cats.data.EitherT.fromEither
import cats.data.{EitherT, NonEmptyList, ValidatedNel}
import cats.effect.IO
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._

import scala.util.Try

object Parse {

  type Parse[A] = EitherT[IO, NonEmptyList[String], A]

  type ParseValidated[A] = IO[ValidatedNel[String, A]]

  def error[A](error: String, tail: String*): Parse[A] = EitherT.leftT {
    NonEmptyList.of(error, tail: _*)
  }


  // If anyone has a magic import that does exactly what these helpers do, please replace, thx!

  def errorOrParse[A](f: => ErrorOr[A]): Parse[A] = fromEither[IO](f.toEither)

  def tryParse[A](f: => Try[A]): Parse[A] = errorOrParse(f.toErrorOr)

  def goParse[A](f: => A): Parse[A] = tryParse(Try(f))

}
