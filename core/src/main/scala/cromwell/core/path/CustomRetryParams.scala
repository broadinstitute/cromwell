package cromwell.core.path

import common.util.Backoff
import cromwell.core.retry.SimpleExponentialBackoff

import scala.concurrent.duration.Duration
import scala.concurrent.duration._
import scala.language.postfixOps

object CustomRetryParams {
  val Default = CustomRetryParams(
    timeout = Duration.Inf,
    maxRetries = Option(3),
    backoff = SimpleExponentialBackoff(1 seconds, 3 seconds, 1.5d),
    isTransient = throwableToFalse,
    isFatal = throwableToFalse
  )

  def throwableToFalse(t: Throwable) = false
}

case class CustomRetryParams(timeout: Duration,
                             maxRetries: Option[Int],
                             backoff: Backoff,
                             isTransient: Throwable => Boolean,
                             isFatal: Throwable => Boolean
)
