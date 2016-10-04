package cromwell.core.path

import cromwell.core.retry.SimpleExponentialBackoff

import scala.concurrent.duration.Duration
import scala.concurrent.duration._
import scala.language.postfixOps

object RetryParams {
  val Default = RetryParams(
    timeout = 1 second,
    maxRetries = Option(3),
    backoff = SimpleExponentialBackoff(1 seconds, 3 seconds, 1.5D),
    isTransient = throwableToFalse,
    isFatal = throwableToFalse
  )

  def throwableToFalse(t: Throwable) = false
}

case class RetryParams(timeout: Duration,
                       maxRetries: Option[Int],
                       backoff: SimpleExponentialBackoff,
                       isTransient: Throwable => Boolean,
                       isFatal: Throwable => Boolean)
