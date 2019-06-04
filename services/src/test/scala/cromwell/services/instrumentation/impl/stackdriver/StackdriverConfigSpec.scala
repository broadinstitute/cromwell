package cromwell.services.instrumentation.impl.stackdriver

import com.typesafe.config.ConfigFactory
import common.exception.AggregatedMessageException
import cromwell.core.TestKitSuite
import org.scalatest.{BeforeAndAfterAll, FlatSpecLike, Matchers}

import scala.concurrent.duration._

class StackdriverConfigSpec extends TestKitSuite with FlatSpecLike with BeforeAndAfterAll with Matchers {
  behavior of "StackdriverConfig"

  val googleConfig = ConfigFactory.parseString(
    s"""
       |google {
       |  application-name = "cromwell"
       |  auths = [
       |    {
       |      name = "application-default"
       |      scheme = "application_default"
       |    }
       |  ]
       |}
      """.stripMargin
  )

  it should "correctly parse all config" in {
    val config = ConfigFactory.parseString(
      """
        |auth = "application-default"
        |google-project = "my-project"
        |flush-rate = 1 minute
        |cromwell-instance-role = "backend"
        |cromwell-perf-test-case = "perf-test-1"
      """.stripMargin
    )

    val stackdriverConfig = StackdriverConfig(config, googleConfig)

    stackdriverConfig.auth.name shouldBe "application-default"
    stackdriverConfig.googleProject shouldBe "my-project"
    stackdriverConfig.flushRate shouldBe 1.minute
    stackdriverConfig.cromwellInstanceRole shouldBe Option("backend")
    stackdriverConfig.cromwellPerfTestCase shouldBe Option("perf-test-1")
  }


  it should "correctly parse config with optional values" in {
    val config = ConfigFactory.parseString(
      """
        |auth = "application-default"
        |google-project = "my-project"
        |flush-rate = 1 minute
      """.stripMargin
    )

    val stackdriverConfig = StackdriverConfig(config, googleConfig)

    stackdriverConfig.auth.name shouldBe "application-default"
    stackdriverConfig.googleProject shouldBe "my-project"
    stackdriverConfig.flushRate shouldBe 1.minute
    stackdriverConfig.cromwellInstanceIdentifier shouldBe None
    stackdriverConfig.cromwellInstanceRole shouldBe None
    stackdriverConfig.cromwellPerfTestCase shouldBe None
  }


  it should "correctly add cromwell instance identifier to config" in {
    val globalConfig = ConfigFactory.parseString(
      s"""
         |google {
         |  application-name = "cromwell"
         |  auths = [
         |    {
         |      name = "application-default"
         |      scheme = "application_default"
         |    }
         |  ]
         |}
         |
         |system.cromwell_id = "cromwell-101"
      """.stripMargin
    )

    val config = ConfigFactory.parseString(
      """
        |auth = "application-default"
        |google-project = "my-project"
        |flush-rate = 1 minute
      """.stripMargin
    )

    val stackdriverConfig = StackdriverConfig(config, globalConfig)

    stackdriverConfig.auth.name shouldBe "application-default"
    stackdriverConfig.googleProject shouldBe "my-project"
    stackdriverConfig.flushRate shouldBe 1.minute
    stackdriverConfig.cromwellInstanceIdentifier shouldBe Option("cromwell-101")
    stackdriverConfig.cromwellInstanceRole shouldBe None
    stackdriverConfig.cromwellPerfTestCase shouldBe None
  }


  it should "throw error for invalid auth" in {
    val config = ConfigFactory.parseString(
      """
        |auth = "my-auth"
        |google-project = "my-project"
        |flush-rate = 1 minute
      """.stripMargin
    )

    val exception = the[AggregatedMessageException] thrownBy StackdriverConfig(config, googleConfig)

    exception.getMessage shouldBe """Stackdriver instrumentation config is invalid. Error(s):
                                    |`auth` scheme is invalid. Errors: NonEmptyList(`google` configuration stanza does not contain an auth named 'my-auth'.  Known auth names: application-default)""".stripMargin
  }


  it should "throw error for invalid flush rate" in {
    val config = ConfigFactory.parseString(
      """
        |auth = "application-default"
        |google-project = "my-project"
        |flush-rate = 30 seconds
      """.stripMargin
    )

    val exception = the[AggregatedMessageException] thrownBy StackdriverConfig(config, googleConfig)

    exception.getMessage shouldBe """Stackdriver instrumentation config is invalid. Error(s):
                                    |`flush-rate` must be 1 minute or longer, specified rate is `30 seconds`. Google's Stackdriver API needs each metric to be sent not more than once per minute.""".stripMargin
  }
}
