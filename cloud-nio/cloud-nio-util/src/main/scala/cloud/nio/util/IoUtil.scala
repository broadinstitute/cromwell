package cloud.nio.util

import cats.data.NonEmptyList
import cats.effect.IO

object IoUtil {

  implicit class EnhancedNonEmptyList[A](val tries: NonEmptyList[IO[A]]) extends AnyVal {

    def trySequence(failureMessage: String): IO[A] = {

      type ACC = IO[Either[NonEmptyList[Exception], A]]

      def attemptHead(headIo: IO[A]): ACC = {
        attemptIo(NonEmptyList.one)(headIo)
      }

      def attemptAcc(accIo: ACC, nextIo: IO[A]): ACC = {
        accIo flatMap {
          case Right(previousSuccess)   => IO.pure(Right(previousSuccess))
          case Left(previousExceptions) => attemptIo(_ :: previousExceptions)(nextIo)
        }
      }

      def attemptIo(f: Exception => NonEmptyList[Exception])(io: IO[A]): ACC = {
        io.attempt flatMap {
          case Right(success)             => IO.pure(Right(success))
          case Left(exception: Exception) => IO.pure(Left(f(exception)))
          case Left(throwable)            => throw throwable
        }
      }

      val res: ACC = tries.tail.foldLeft(attemptHead(tries.head))(attemptAcc)

      val eitherIo: IO[Either[Throwable, A]] = res map {
        _.left.map(nel => new AggregateException(failureMessage, nel.reverse))
      }

      eitherIo.flatMap(IO.fromEither)
    }

  }
}
