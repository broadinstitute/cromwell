package centaur.api

import akka.actor.ActorSystem
import akka.pattern.after

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}

object Retry {

  /**
    * Copied from cromwell.core
    * Replaced the backoff with a fixed retry delay
    */
  def withRetry[A](f: () => Future[A],
                   maxRetries: Option[Int] = Option(10),
                   delay: FiniteDuration,
                   isTransient: Throwable => Boolean = throwableToFalse,
                   isFatal: Throwable => Boolean = throwableToFalse,
                   onRetry: Throwable => Unit = noopOnRetry
  )(implicit actorSystem: ActorSystem): Future[A] = {
    // In the future we might want EC passed in separately but at the moment it caused more issues than it solved to do so
    implicit val ec: ExecutionContext = actorSystem.dispatcher

    f() recoverWith {
      case throwable if isFatal(throwable) => Future.failed(throwable)
      case throwable if !isFatal(throwable) =>
        val retriesLeft = if (isTransient(throwable)) maxRetries else maxRetries map { _ - 1 }

        if (retriesLeft.forall(_ > 0)) {
          onRetry(throwable)
          after(delay, actorSystem.scheduler)(
            withRetry(f,
                      delay = delay,
                      maxRetries = retriesLeft,
                      isTransient = isTransient,
                      isFatal = isFatal,
                      onRetry = onRetry
            )
          )
        } else {
          Future.failed(throwable)
        }
    }
  }

  def throwableToFalse(t: Throwable) = false
  def noopOnRetry(t: Throwable) = {}
}
