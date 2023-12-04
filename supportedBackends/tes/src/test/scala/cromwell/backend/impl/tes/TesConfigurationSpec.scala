package cromwell.backend.impl.tes

import com.typesafe.config.{Config, ConfigException, ConfigFactory}
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.core.retry.SimpleExponentialBackoff
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.concurrent.duration._
import scala.language.postfixOps

class TesConfigurationSpec extends AnyFlatSpec with Matchers {

  behavior of "TesConfiguration"

  def makeTesConfig(backendConfig: Config) = new TesConfiguration(
    BackendConfigurationDescriptor(
      backendConfig,
      ConfigFactory.empty()
    )
  )

  def backoffsAreEquivalent(expectedBackoff: SimpleExponentialBackoff,
                            actualBackoff: SimpleExponentialBackoff
  ): Boolean = {
    val b1 = expectedBackoff.googleBackoff
    val b2 = actualBackoff.googleBackoff
    b1.getInitialIntervalMillis == b2.getInitialIntervalMillis &&
    b1.getMaxIntervalMillis == b2.getMaxIntervalMillis &&
    b1.getMultiplier == b2.getMultiplier &&
    b1.getRandomizationFactor == b2.getRandomizationFactor
  }

  it should "use default backoffs when no custom config provided" in {
    val tesConfig = makeTesConfig(TesTestConfig.backendConfig)
    backoffsAreEquivalent(TesConfiguration.defaultPollBackoff, tesConfig.pollBackoff) shouldBe true
    backoffsAreEquivalent(TesConfiguration.defaultExecOrRecoverBackoff, tesConfig.executeOrRecoverBackoff) shouldBe true
  }

  it should "use configured backoffs if they exist" in {
    val tesConfig = makeTesConfig(TesTestConfig.backendConfigWithBackoffs)
    backoffsAreEquivalent(SimpleExponentialBackoff(5 seconds, 1 minute, 2.5, .7), tesConfig.pollBackoff) shouldBe true
    backoffsAreEquivalent(SimpleExponentialBackoff(3 minutes, 1 hours, 5, .1),
                          tesConfig.executeOrRecoverBackoff
    ) shouldBe true
  }

  it should "fail if user defines an invalid backoff" in {
    assertThrows[ConfigException.Missing](makeTesConfig(TesTestConfig.backendConfigWithInvalidBackoffs))
  }
}
