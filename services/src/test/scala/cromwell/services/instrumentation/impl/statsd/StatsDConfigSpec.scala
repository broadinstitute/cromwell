package cromwell.services.instrumentation.impl.statsd

import com.typesafe.config.ConfigFactory
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.concurrent.duration._

class StatsDConfigSpec extends AnyFlatSpec with Matchers {
  behavior of "StatsDConfig"
  
  it should "parse correct service configuration" in {
    val config = ConfigFactory.parseString(
      """
        |hostname = "localhost"
        |port = 8125
        |prefix = "prefix_value" # can be used to prefix all metrics with an api key for example
        |flush-rate = 1 second # rate at which aggregated metrics will be sent to statsd
      """.stripMargin
    )

    val statsDConfig = StatsDConfig(config)

    statsDConfig.hostname shouldBe "localhost"
    statsDConfig.port shouldBe 8125
    statsDConfig.prefix shouldBe Some("prefix_value")
    statsDConfig.flushRate shouldBe 1.second
  }

  it should "not parse incorrect service configuration" in {
    val config = ConfigFactory.parseString("{}")

    val exception = the[IllegalArgumentException] thrownBy StatsDConfig(config)

    exception.getMessage shouldBe """StatsD config is invalid:
                                    |String: 1: No configuration setting found for key 'hostname'
                                    |String: 1: No configuration setting found for key 'port'
                                    |String: 1: No configuration setting found for key 'flush-rate'""".stripMargin
  }

}
