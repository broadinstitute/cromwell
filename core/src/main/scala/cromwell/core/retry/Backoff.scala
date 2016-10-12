package cromwell.core.retry

import com.google.api.client.util.ExponentialBackOff

import scala.concurrent.duration.{Duration, FiniteDuration}

trait Backoff {
  /** Next interval in millis */
  def backoffMillis: Long
  /** Get the next instance of backoff. This should be called after every call to backoffMillis */
  def next: Backoff
}

object InitialGapBackoff {
  def apply(initialGap: FiniteDuration, initialInterval: FiniteDuration, maxInterval: FiniteDuration, multiplier: Double) = {
    new InitialGapBackoff(initialGap, new ExponentialBackOff.Builder()
      .setInitialIntervalMillis(initialInterval.toMillis.toInt)
      .setMaxIntervalMillis(maxInterval.toMillis.toInt)
      .setMultiplier(multiplier)
      .setMaxElapsedTimeMillis(Int.MaxValue)
      .build())
  }
}

case class InitialGapBackoff(initialGapMillis: FiniteDuration, googleBackoff: ExponentialBackOff) extends Backoff {
  assert(initialGapMillis.compareTo(Duration.Zero) != 0, "Initial gap cannot be null, use SimpleBackoff instead.")

  override val backoffMillis = initialGapMillis.toMillis
  /** Switch to a SimpleExponentialBackoff after the initial gap has been used */
  override def next = new SimpleExponentialBackoff(googleBackoff)
}

object SimpleExponentialBackoff {
  def apply(initialInterval: FiniteDuration, maxInterval: FiniteDuration, multiplier: Double) = {
    new SimpleExponentialBackoff(new ExponentialBackOff.Builder()
      .setInitialIntervalMillis(initialInterval.toMillis.toInt)
      .setMaxIntervalMillis(maxInterval.toMillis.toInt)
      .setMultiplier(multiplier)
      .setMaxElapsedTimeMillis(Int.MaxValue)
      .build())
  }
}

case class SimpleExponentialBackoff(googleBackoff: ExponentialBackOff) extends Backoff {
  override def backoffMillis = googleBackoff.nextBackOffMillis()
  /** google ExponentialBackOff is mutable so we can keep returning the same instance */
  override def next = this
}
