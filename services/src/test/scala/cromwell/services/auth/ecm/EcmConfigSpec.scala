package cromwell.services.auth.ecm

import com.typesafe.config.ConfigFactory
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class EcmConfigSpec extends AnyFlatSpec with Matchers {

  it should "parse ECM base url when present" in {
    val config = ConfigFactory.parseString(s"""
                                              |enabled = true
                                              |auth.azure = true
                                              |ecm.base-url = "https://mock-ecm-url.org"
      """.stripMargin)

    val actualEcmConfig = EcmConfig(config)

    actualEcmConfig.baseUrl shouldBe defined
    actualEcmConfig.baseUrl.get shouldBe "https://mock-ecm-url.org"
  }

  it should "return None when ECM base url is absent" in {
    val config = ConfigFactory.parseString(s"""
                                              |enabled = true
                                              |auth.azure = true
      """.stripMargin)

    EcmConfig(config).baseUrl shouldBe None
  }
}
