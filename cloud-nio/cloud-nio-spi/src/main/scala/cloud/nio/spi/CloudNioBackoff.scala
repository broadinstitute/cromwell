package cloud.nio.spi

import com.google.api.client.util.ExponentialBackOff

import scala.concurrent.duration.{Duration, FiniteDuration}

trait CloudNioBackoff {
  /** Next interval in millis */
  def backoffMillis: Long
  /** Get the next instance of backoff. This should be called after every call to backoffMillis */
  def next: CloudNioBackoff
}

object CloudNioInitialGapBackoff {
  def apply(initialGap: FiniteDuration, initialInterval: FiniteDuration, maxInterval: FiniteDuration, multiplier: Double) = {
    new CloudNioInitialGapBackoff(initialGap, new ExponentialBackOff.Builder()
      .setInitialIntervalMillis(initialInterval.toMillis.toInt)
      .setMaxIntervalMillis(maxInterval.toMillis.toInt)
      .setMultiplier(multiplier)
      .setMaxElapsedTimeMillis(Int.MaxValue)
      .build())
  }
}

case class CloudNioInitialGapBackoff(initialGapMillis: FiniteDuration, googleBackoff: ExponentialBackOff) extends CloudNioBackoff {
  assert(initialGapMillis.compareTo(Duration.Zero) != 0, "Initial gap cannot be null, use SimpleBackoff instead.")

  override val backoffMillis = initialGapMillis.toMillis
  /** Switch to a SimpleExponentialBackoff after the initial gap has been used */
  override def next = new CloudNioSimpleExponentialBackoff(googleBackoff)
}

object CloudNioSimpleExponentialBackoff {
  def apply(initialInterval: FiniteDuration, maxInterval: FiniteDuration, multiplier: Double) = {
    new CloudNioSimpleExponentialBackoff(new ExponentialBackOff.Builder()
      .setInitialIntervalMillis(initialInterval.toMillis.toInt)
      .setMaxIntervalMillis(maxInterval.toMillis.toInt)
      .setMultiplier(multiplier)
      .setMaxElapsedTimeMillis(Int.MaxValue)
      .build())
  }
}

case class CloudNioSimpleExponentialBackoff(googleBackoff: ExponentialBackOff) extends CloudNioBackoff {
  override def backoffMillis = googleBackoff.nextBackOffMillis()
  /** google ExponentialBackOff is mutable so we can keep returning the same instance */
  override def next = this
}
