package cloud.nio.spi

import cloud.nio.spi.CloudNioBackoff._
import com.google.api.client.util.ExponentialBackOff

import scala.concurrent.duration.{Duration, FiniteDuration}

trait CloudNioBackoff {
  /** Next interval in millis */
  def backoffMillis: Long
  /** Get the next instance of backoff. This should be called after every call to backoffMillis */
  def next: CloudNioBackoff
}

object CloudNioBackoff {
  private[spi] def newExponentialBackOff(initialInterval: FiniteDuration,
                                         maxInterval: FiniteDuration,
                                         multiplier: Double,
                                         randomizationFactor: Double,
                                        ): ExponentialBackOff = {
    new ExponentialBackOff.Builder()
      .setInitialIntervalMillis(initialInterval.toMillis.toInt)
      .setMaxIntervalMillis(maxInterval.toMillis.toInt)
      .setMultiplier(multiplier)
      .setRandomizationFactor(randomizationFactor)
      .setMaxElapsedTimeMillis(Int.MaxValue)
      .build()
  }
}

object CloudNioInitialGapBackoff {
  def apply(initialGap: FiniteDuration,
            initialInterval: FiniteDuration,
            maxInterval: FiniteDuration,
            multiplier: Double,
            randomizationFactor: Double = ExponentialBackOff.DEFAULT_RANDOMIZATION_FACTOR,
           ): CloudNioInitialGapBackoff = {
    new CloudNioInitialGapBackoff(
      initialGap,
      newExponentialBackOff(
        initialInterval = initialInterval,
        maxInterval = maxInterval,
        multiplier = multiplier,
        randomizationFactor = randomizationFactor,
      )
    )
  }
}

case class CloudNioInitialGapBackoff(initialGapMillis: FiniteDuration, googleBackoff: ExponentialBackOff) extends CloudNioBackoff {
  assert(initialGapMillis.compareTo(Duration.Zero) != 0, "Initial gap cannot be null, use SimpleBackoff instead.")

  override val backoffMillis: Long = initialGapMillis.toMillis
  /** Switch to a SimpleExponentialBackoff after the initial gap has been used */
  override def next = new CloudNioSimpleExponentialBackoff(googleBackoff)
}

object CloudNioSimpleExponentialBackoff {
  def apply(initialInterval: FiniteDuration,
            maxInterval: FiniteDuration,
            multiplier: Double,
            randomizationFactor: Double = ExponentialBackOff.DEFAULT_RANDOMIZATION_FACTOR,
           ): CloudNioSimpleExponentialBackoff = {
    new CloudNioSimpleExponentialBackoff(
      newExponentialBackOff(
        initialInterval = initialInterval,
        maxInterval = maxInterval,
        multiplier = multiplier,
        randomizationFactor = randomizationFactor,
      )
    )
  }
}

case class CloudNioSimpleExponentialBackoff(googleBackoff: ExponentialBackOff) extends CloudNioBackoff {
  override def backoffMillis: Long = googleBackoff.nextBackOffMillis()
  /** google ExponentialBackOff is mutable so we can keep returning the same instance */
  override def next: CloudNioBackoff = this
}
