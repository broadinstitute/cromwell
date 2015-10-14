package lenthall.config

import java.net.URL

import com.typesafe.config.ConfigFactory
import lenthall.config.ValidatedConfig._
import lenthall.test.logging.TestLogger
import org.scalatest.{FlatSpec, Matchers}

import scalaz.Scalaz._

class ValidatedConfigSpec extends FlatSpec with Matchers {
  behavior of "ValidatedConfig"

  val exampleConfig = ConfigFactory.parseString(
    """
      |{
      |  stringVal: "string"
      |  booleanVal: true
      |  intVal: 123
      |  longVal: 123
      |  doubleVal: 12.3
      |  urlVal: "https://example.org"
      |  urlBad: "htts://example.org"
      |}
    """.stripMargin)

  val missingNel = "Could not find key: missing".failureNel

  it should "return the config value as a success nel when present" in {
    exampleConfig.validateString("stringVal") should be("string".successNel)
    exampleConfig.validateBoolean("booleanVal") should be(true.successNel)
    exampleConfig.validateInt("intVal") should be(123.successNel)
    exampleConfig.validateLong("longVal") should be(123L.successNel)
    exampleConfig.validateDouble("doubleVal") should be(12.3.successNel)
  }

  it should "return a failure nel when missing" in {
    exampleConfig.validateString("missing") should be(missingNel)
    exampleConfig.validateBoolean("missing") should be(missingNel)
    exampleConfig.validateInt("missing") should be(missingNel)
    exampleConfig.validateLong("missing") should be(missingNel)
    exampleConfig.validateDouble("missing") should be(missingNel)
  }

  it should "check urls correctly" in {
    exampleConfig.validateURL("urlVal") should be(new URL("https://example.org").successNel)
    exampleConfig.validateURL("urlBad") should be("unknown protocol: htts".failureNel)
    exampleConfig.validateURL("missing") should be(missingNel)
  }

  it should "log the unexpected config keys" in {
    import TestLogger._
    import ValidatedConfig._
    val config = ConfigFactory.parseString("{ unexpectedKey1 = true, unexpectedKey2 = 1 }")
    withTestLoggerFor(validationLogger) { logger =>
      config.warnNotRecognized(Set("expectedKey"), "testContext")
      logger.messages should be(
        """|[WARN] Unrecognized configuration key(s) for testContext: unexpectedKey1, unexpectedKey2
           |""".stripMargin)
    }
  }

  it should "not log the expected config keys" in {
    import TestLogger._
    import ValidatedConfig._
    val config = ConfigFactory.parseString("{ expectedKey1 = true, unexpectedKey2 = 1 }")
    withTestLoggerFor(validationLogger) { logger =>
      config.warnNotRecognized(Set("expectedKey1"), "testContext")
      logger.messages should be(
        """|[WARN] Unrecognized configuration key(s) for testContext: unexpectedKey2
           |""".stripMargin)
    }
  }

  it should "not log any keys when all match" in {
    import TestLogger._
    import ValidatedConfig._
    val config = ConfigFactory.parseString("{ expectedKey1 = true, expectedKey2 = 1 }")
    withTestLoggerFor(validationLogger) { logger =>
      config.warnNotRecognized(Set("expectedKey1", "expectedKey2"), "testContext")
      logger.messages should be("")
    }
  }
}
