package cromwell.util

import java.sql.SQLTransientException

import akka.actor.ActorSystem
import cromwell.core.retry.{Retry, SimpleExponentialBackoff}

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.language.postfixOps

object DatabaseUtil {
  private def isTransient(throwable: Throwable): Boolean = throwable match {
    case _: SQLTransientException => true
    case _ => false
  }

  def withRetry[A](f: () => Future[A])(implicit actorSystem: ActorSystem): Future[A] = {
    val RetryBackoff = SimpleExponentialBackoff(50 millis, 1 seconds, 1D)
    Retry.withRetry(f, maxRetries = Option(10), backoff = RetryBackoff, isTransient = isTransient)
  }
}
