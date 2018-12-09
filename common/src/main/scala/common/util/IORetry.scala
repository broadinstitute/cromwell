package common.util

import cats.effect.{IO, Timer}

import scala.concurrent.duration._
import scala.util.control.NonFatal

object IORetry {
  def noOpOnRetry[S]: (Throwable, S) => S = (_, s) => s
  
  object StatefulIoError {
    def noop[S] = new StatefulIoError[S] {
      override def toThrowable(state: S, throwable: Throwable) = throwable
    }
  }

  implicit val noopUnitState = IORetry.StatefulIoError.noop[Unit]

  /**
    * When we reach a point where we need to fail the IO (because we ran out of retries, or exception was fatal etc...)
    * This provides a method to generate a Throwable from the current state and the last failure.
    */
  trait StatefulIoError[S] {
    def toThrowable(state: S, throwable: Throwable): Throwable
  }

  /**
    * Same as Retry.withRetry using cats.effect.IO
    * The state is being passed down upon each retry. It can be use to accumulate information as IO operations
    * fail and are being retried.
    */
  def withRetry[A, S](io: IO[A],
                      state: S,
                      maxRetries: Option[Int],
                      backoff: Backoff,
                      isTransient: Throwable => Boolean = throwableToFalse,
                      isFatal: Throwable => Boolean = throwableToFalse,
                      onRetry: (Throwable, S) => S = noOpOnRetry[S])
                     (implicit timer: Timer[IO], statefulIoException: StatefulIoError[S]): IO[A] = {
    lazy val delay = backoff.backoffMillis.millis

    def fail(throwable: Throwable) = IO.raiseError(statefulIoException.toThrowable(state, throwable))

    io handleErrorWith {
      case throwable if isFatal(throwable) => fail(throwable)
      case NonFatal(throwable) =>
        val retriesLeft = if (isTransient(throwable)) maxRetries else maxRetries map { _ - 1 }

        if (retriesLeft.forall(_ > 0)) {
          for {
            _ <- IO.sleep(delay)
            retried <- withRetry(io, onRetry(throwable, state), retriesLeft, backoff.next, isTransient, isFatal, onRetry)
          } yield retried
        }
        else fail(throwable)
      case fatal => throw fatal
    }

  }

  def throwableToFalse(t: Throwable) = false
}
