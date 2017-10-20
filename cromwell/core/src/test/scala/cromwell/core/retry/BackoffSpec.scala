package cromwell.core.retry

import com.google.api.client.util.ExponentialBackOff.Builder
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.duration._

class BackoffSpec extends FlatSpec with Matchers {

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

}
