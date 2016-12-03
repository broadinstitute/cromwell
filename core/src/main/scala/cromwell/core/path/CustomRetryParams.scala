package cromwell.core.path

import cromwell.core.retry.{Backoff, SimpleExponentialBackoff}

import scala.concurrent.duration.Duration
import scala.concurrent.duration._
import scala.language.postfixOps

object CustomRetryParams {
  val Default = CustomRetryParams(
    timeout = Duration.Inf,
    maxRetries = Option(3),
    backoff = SimpleExponentialBackoff(1 seconds, 3 seconds, 1.5D),
    isTransient = throwableToFalse,
    isFatal = throwableToFalse
  )

  def throwableToFalse(t: Throwable) = false
}

case class CustomRetryParams(timeout: Duration,
                             maxRetries: Option[Int],
                             backoff: Backoff,
                             isTransient: Throwable => Boolean,
                             isFatal: Throwable => Boolean)
