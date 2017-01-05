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
    * @param f A function Unit => Future which will be executed once per cycle.
    * @param maxRetries An optional number of times to retry the thunk. If this is None, there is no limit.
    * @param backoff An exponential backoff strategy to use for each retry strategy.
    * @param isTransient An optional function of Throwable => Boolean. If provided and it returns true, this throwable
    *                    does not count towards the maxRetries limit.
    * @param isFatal An optional function of Throwable => Boolean. If provided and it returns true, withRetry will halt
    * @param actorSystem An ActorSystem, used for its scheduler and dispatcher.
    * @tparam A The return type of the thunk
    * @return The final completed Future
    */
  def withRetry[A](f: () => Future[A],
                   maxRetries: Option[Int] = Option(10),
                   backoff: Backoff = SimpleExponentialBackoff(5 seconds, 10 seconds, 1.1D),
                   isTransient: Throwable => Boolean = throwableToFalse,
                   isFatal: Throwable => Boolean = throwableToFalse)
                   (implicit actorSystem: ActorSystem): Future[A] = {
    // In the future we might want EC passed in separately but at the moment it caused more issues than it solved to do so
    implicit val ec: ExecutionContext = actorSystem.dispatcher
    val delay = backoff.backoffMillis.millis

    if (maxRetries.forall(_ > 0)) {
      f() recoverWith {
        case throwable if isFatal(throwable) => Future.failed(CromwellFatalException(throwable))
        case throwable if !isFatal(throwable) =>
          val retriesLeft = if (isTransient(throwable)) maxRetries else maxRetries map { _ - 1 }
          after(delay, actorSystem.scheduler)(withRetry(f, backoff = backoff, maxRetries = retriesLeft, isTransient = isTransient, isFatal = isFatal))
      }
    } else f() recoverWith {
      case e: Exception => Future.failed(new CromwellFatalException(e))
    }
  }

  def throwableToFalse(t: Throwable) = false
}

