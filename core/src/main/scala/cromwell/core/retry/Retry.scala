package cromwell.core.retry

import akka.actor.ActorSystem
import cromwell.core.CromwellFatalException

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.concurrent.{ExecutionContext, Future}
import akka.pattern.after

object Retry {
  /**
    * Retries a Future on a designated backoff strategy until either a designated number of retries or a fatal error
    * is reached.
    *
    * @param f A Future thunk which will be executed once per cycle
    * @param maxRetries An optional number of times to retry the thunk. If this is None, there is no limit.
    * @param backoff An exponential backoff strategy to use for each retry strategy.
    * @param isTransient An optional function of Throwable => Boolean. If provided and it returns true, this throwable
    *                    does not count towards the maxRetries limit.
    * @param isFatal An optional function of Throwable => Boolean. If provided and it returns true, withRetry will halt
    * @param actorSystem An ActorSystem, used for its scheduler and dispatcher.
    * @tparam A The return type of the thunk
    * @return The final completed Future
    */
  def withRetry[A](f: => Future[A],
                   maxRetries: Option[Int] = Option(10),
                   backoff: SimpleExponentialBackoff = SimpleExponentialBackoff(5 seconds, 10 seconds, 1.1D),
                   isTransient: Option[Throwable => Boolean] = None,
                   isFatal: Option[Throwable => Boolean] = None)
                   (implicit actorSystem: ActorSystem): Future[A] = {
    // In the future we might want EC passed in separately but at the moment it caused more issues than it solved to do so
    implicit val ec: ExecutionContext = actorSystem.dispatcher
    val delay = backoff.backoffMillis.millis

    if (maxRetries.forall(_ > 0)) {
      f recoverWith {
        case throwable if isFatal.evaluate(throwable) => Future.failed(new CromwellFatalException(throwable))
        case throwable if !isFatal.evaluate(throwable) =>
          val retriesLeft = if (isTransient.evaluate(throwable)) maxRetries else maxRetries map { _ - 1 }
          after(delay, actorSystem.scheduler)(withRetry(f, backoff = backoff, maxRetries = retriesLeft))
      }
    } else f recoverWith {
      case e: Exception => Future.failed(new CromwellFatalException(e))
    }
  }

  implicit class EnhancedOptionThrowableToBoolean(val f: Option[Throwable => Boolean]) extends AnyVal {
    /** If f is defined, run the function using the passed throwable otherwise return false */
    def evaluate(throwable: Throwable) = f.fold(false) { _(throwable) }
  }
}

