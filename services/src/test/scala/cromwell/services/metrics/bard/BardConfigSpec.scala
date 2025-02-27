package cromwell.services.metrics.bard

import com.typesafe.config.ConfigFactory
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class BardConfigSpec extends AnyFlatSpec with Matchers {

  behavior of "BardConfig"

  it should "parse the config" in {
    val config = ConfigFactory.parseString(s"""
                                              |enabled = true
                                              |bard.base-url = "https://mock-bard-url.org"
                                              |bard.connection-pool-size = 10
      """.stripMargin)

    val actualBardConfig = BardConfig(config)

    actualBardConfig.enabled shouldBe true
    actualBardConfig.baseUrl shouldBe "https://mock-bard-url.org"
    actualBardConfig.connectionPoolSize shouldBe 10
  }

  it should "not enable Bard if its disabled in the config" in {
    val config = ConfigFactory.parseString(s"""
                                              |enabled = false
                                              |bard.base-url = ""
                                              |bard.connection-pool-size = 0
      """.stripMargin)

    val actualBardConfig = BardConfig(config)

    actualBardConfig.enabled shouldBe false
  }
}
