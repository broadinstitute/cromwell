package cromwell.engine.backend.jes

import java.net.URL

import com.typesafe.config.ConfigFactory
import org.scalatest.{Matchers, FlatSpec}
import wdl4s.ThrowableWithErrors

class JesAttributesSpec extends FlatSpec with Matchers {

  it should "parse correct JES config" in {
    val configString = """
          backend {
           jes {
             project = "myProject"
             baseExecutionBucket = "gs://myBucket"
             endpointUrl = "http://myEndpoint"
             maximumPollingInterval = 600
             [PREEMPTIBLE]
           }
          }""".stripMargin

    val fullConfig = ConfigFactory.parseString(configString.replace("[PREEMPTIBLE]", "preemptible = 3"))

    val fullAttributes = JesAttributes.apply(fullConfig)
    fullAttributes.endpointUrl.toString shouldBe new URL("http://myEndpoint").toString
    fullAttributes.project shouldBe "myProject"
    fullAttributes.executionBucket shouldBe "gs://myBucket"
    fullAttributes.maxPollingInterval shouldBe 600

    val noPreemptibleConfig = ConfigFactory.parseString(configString.replace("[PREEMPTIBLE]", ""))

    val noPreemptibleAttributes = JesAttributes.apply(noPreemptibleConfig)
    noPreemptibleAttributes.endpointUrl.toString shouldBe new URL("http://myEndpoint").toString
    noPreemptibleAttributes.project shouldBe "myProject"
    noPreemptibleAttributes.executionBucket shouldBe "gs://myBucket"
    noPreemptibleAttributes.maxPollingInterval shouldBe 600
  }

  it should "not parse invalid config" in {
    val nakedConfig =
      ConfigFactory.parseString("""
        |backend {
        | jes {
        |   endpointUrl = "myEndpoint"
        | }
        |}
      """.stripMargin)

    val exception = intercept[IllegalArgumentException with ThrowableWithErrors] {
      JesAttributes.apply(nakedConfig)
    }
    val errorsList = exception.errors.list
    errorsList should contain ("Could not find key: project")
    errorsList should contain ("Could not find key: baseExecutionBucket")
    errorsList should contain ("no protocol: myEndpoint")
  }

}
