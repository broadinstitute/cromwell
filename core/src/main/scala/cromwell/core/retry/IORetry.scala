package cromwell.core.retry

import cats.effect.{IO, Timer}
import cats.syntax.all._

import scala.concurrent.duration._

object IORetry {
  def noOpOnRetry[S]: (Throwable, S) => S = (_, s) => s

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
                      backoff: Backoff = SimpleExponentialBackoff(5.seconds, 10.seconds, 1.1D),
                      isTransient: Throwable => Boolean = throwableToFalse,
                      isFatal: Throwable => Boolean = throwableToFalse,
                      onRetry: (Throwable, S) => S = noOpOnRetry)
                     (implicit timer: Timer[IO], statefulIoException: StatefulIoError[S]): IO[A] = {
    val delay = backoff.backoffMillis.millis

    def fail(throwable: Throwable) = IO.raiseError(statefulIoException.toThrowable(state, throwable))

    io handleErrorWith {
      case throwable if isFatal(throwable) => fail(throwable)
      case throwable =>
        val retriesLeft = if (isTransient(throwable)) maxRetries else maxRetries map { _ - 1 }

        if (retriesLeft.forall(_ > 0)) {
          IO.sleep(delay) *> withRetry(io, onRetry(throwable, state), retriesLeft, backoff.next, isTransient, isFatal, onRetry)
        }
        else fail(throwable)
    }

  }

  def throwableToFalse(t: Throwable) = false
}
