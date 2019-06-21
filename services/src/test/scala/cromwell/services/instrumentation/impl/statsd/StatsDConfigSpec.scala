package cromwell.services.instrumentation.impl.statsd

import com.typesafe.config.ConfigFactory
import org.scalatest.{FlatSpec, Matchers}
import scala.concurrent.duration._

class StatsDConfigSpec extends FlatSpec with Matchers {
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
                                    |No configuration setting found for key 'hostname'
                                    |No configuration setting found for key 'port'
                                    |No configuration setting found for key 'flush-rate'""".stripMargin
  }

}
