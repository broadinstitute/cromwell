package cromwell.core.retry

import java.sql.SQLTransactionRollbackException

import akka.actor.ActorSystem
import akka.pattern.after
import com.typesafe.scalalogging.StrictLogging
import common.util.Backoff
import cromwell.core.CromwellFatalException

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps

object Retry extends StrictLogging {

  def throwableToFalse(t: Throwable) = false

  def noopOnRetry(t: Throwable) = {}

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
                   isFatal: Throwable => Boolean = throwableToFalse,
                   onRetry: Throwable => Unit = noopOnRetry)
                   (implicit actorSystem: ActorSystem): Future[A] = {
    // In the future we might want EC passed in separately but at the moment it caused more issues than it solved to do so
    implicit val ec: ExecutionContext = actorSystem.dispatcher
    val delay = backoff.backoffMillis.millis

    f() recoverWith {
      case throwable if isFatal(throwable) => Future.failed(CromwellFatalException(throwable))
      case throwable if !isFatal(throwable) =>
        val retriesLeft = if (isTransient(throwable)) maxRetries else maxRetries map { _ - 1 }
        
        if (retriesLeft.forall(_ > 0)) {
          onRetry(throwable)
          after(delay, actorSystem.scheduler)(withRetry(f, backoff = backoff.next, maxRetries = retriesLeft, isTransient = isTransient, isFatal = isFatal, onRetry = onRetry))
        } else {
          Future.failed(new CromwellFatalException(throwable))
        }
    }
  }

  /**
   * Retries a Future if a 'SQLTransactionRollbackException' (i.e. deadlock exception) occurs on a designated
   * backoff strategy until a designated number of retries is reached.
   *
   * @param f A function Unit => Future which will be executed once per cycle.
   * @param maxRetries Number of times to retry the future
   * @param backoff An exponential backoff strategy to use for each retry strategy.
   * @tparam A The return type of the Future
   * @return The final completed Future[A]
   */
  def withRetryForTransactionRollback[A](f: () => Future[A],
                                         maxRetries: Int = 5,
                                         backoff: Backoff = SimpleExponentialBackoff(5 seconds, 10 seconds, 1.1D))
                                        (implicit actorSystem: ActorSystem, ec: ExecutionContext): Future[A] = {
    val delay = backoff.backoffMillis.millis

    f() recoverWith {
      case throwable if throwable.isInstanceOf[SQLTransactionRollbackException] =>
        val retriesLeft = maxRetries - 1
        if (retriesLeft > 0) {
          logger.info(s"Received 'SQLTransactionRollbackException'. Retries left $retriesLeft. Will retry now...")
          after(delay, actorSystem.scheduler)(withRetryForTransactionRollback(f, retriesLeft, backoff.next))
        } else {
          Future.failed(new CromwellFatalException(throwable))
        }
    }
  }
}

