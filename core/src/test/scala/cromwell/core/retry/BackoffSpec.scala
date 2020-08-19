package cromwell.core.retry

import com.google.api.client.util.ExponentialBackOff.Builder
import com.typesafe.config.ConfigFactory
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.collection.JavaConverters._
import scala.concurrent.duration._

class BackoffSpec extends AnyFlatSpec with Matchers {

  it should "honor initial gap" in {
    val exponentialBackoff = new InitialGapBackoff(
      3.seconds,
      new Builder()
        .setInitialIntervalMillis(1.second.toMillis.toInt)
        .setMaxIntervalMillis(2.seconds.toMillis.toInt)
        .setMaxElapsedTimeMillis(Integer.MAX_VALUE)
        .setRandomizationFactor(0D)
        .build()
    )


    exponentialBackoff.backoffMillis shouldBe 3.seconds.toMillis
    exponentialBackoff.next.backoffMillis shouldBe 1.second.toMillis
  }

  it should "work for simple backoffs" in {
    new SimpleExponentialBackoff(
      new Builder()
        .setInitialIntervalMillis(1.second.toMillis.toInt)
        .setMaxIntervalMillis(2.seconds.toMillis.toInt)
        .setMaxElapsedTimeMillis(Integer.MAX_VALUE)
        .setRandomizationFactor(0D)
        .build()
    ).backoffMillis shouldBe 1.second.toMillis
  }

  it should "throw an exception when trying create a initial gap equal to 0" in {
    an[AssertionError] should be thrownBy {
      new InitialGapBackoff(
        0.millis,
        new Builder()
          .setInitialIntervalMillis(1.second.toMillis.toInt)
          .setMaxIntervalMillis(2.seconds.toMillis.toInt)
          .setMaxElapsedTimeMillis(Integer.MAX_VALUE)
          .setRandomizationFactor(0D)
          .build()
      )
    }
  }

  it should "parse config" in {
    val config = ConfigFactory.parseMap(
      Map(
        "min" -> "5 seconds",
        "max" -> "30 seconds",
        "multiplier" -> 6D,
        "randomization-factor" -> 0D
      ).asJava
    )
    
    val backoff = SimpleExponentialBackoff(config)
    backoff.googleBackoff.getCurrentIntervalMillis shouldBe 5.seconds.toMillis.toInt
    backoff.googleBackoff.getMaxIntervalMillis shouldBe 30.seconds.toMillis.toInt
    backoff.googleBackoff.getMultiplier shouldBe 6D
    backoff.googleBackoff.getRandomizationFactor shouldBe 0D
  }

}
