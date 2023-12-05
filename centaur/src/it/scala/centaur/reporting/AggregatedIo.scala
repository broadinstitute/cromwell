package centaur.reporting

import cats.data.Validated._
import cats.data.ValidatedNel
import cats.effect.IO
import cats.instances.list._
import cats.syntax.either._
import cats.syntax.traverse._

/**
  * Validation that aggregates multiple throwable errors.
  */
object AggregatedIo {

  /**
    * Similar to common.validation.ErrorOr#ErrorOr, but retains the stack traces.
    */
  type ThrowableNelOr[+A] = ValidatedNel[Throwable, A]

  /**
    * Turn a list of IO "inside out" while retaining multiple exceptions within an AggregatedException.
    *
    * @param multipleExceptionContext For two or more exceptions, wrap the exception with this additional context.
    * @param listIo The original IO
    * @tparam A The type of elements stored in the list.
    * @return The "inside out" list.
    */
  def aggregateExceptions[A](multipleExceptionContext: String, listIo: List[IO[A]]): IO[List[A]] = {
    val ioListThrowables: IO[List[ThrowableNelOr[A]]] =
      listIo.traverse[IO, ThrowableNelOr[A]](_.attempt.map(_.toValidatedNel))
    val ioThrowablesList: IO[ThrowableNelOr[List[A]]] = ioListThrowables.map(_.sequence[ThrowableNelOr, A])
    val ioEither: IO[Either[Throwable, List[A]]] = ioThrowablesList map {
      case Valid(list) => Either.right(list)
      case Invalid(nel) if nel.tail.isEmpty => Either.left(nel.head)
      case Invalid(nel) => Either.left(new AggregatedException(multipleExceptionContext, nel.toList))
    }
    ioEither.flatMap(IO.fromEither)
  }

  /**
    * Creates an aggregated exception for multiple exceptions.
    */
  class AggregatedException private[reporting] (exceptionContext: String, suppressed: List[Throwable])
      extends RuntimeException(
        {
          val suppressedZipped = suppressed.zipWithIndex
          val messages = suppressedZipped map { case (throwable, index) =>
            s"\n  ${index + 1}: ${throwable.getMessage}"
          }
          s"$exceptionContext:$messages"
        }
      ) {
    suppressed foreach addSuppressed
  }

}
