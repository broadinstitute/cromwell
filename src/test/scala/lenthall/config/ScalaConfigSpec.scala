package lenthall.config

import com.typesafe.config.ConfigFactory
import lenthall.config.ScalaConfig._
import org.scalatest.{FlatSpec, Matchers}

class ScalaConfigSpec extends FlatSpec with Matchers {
  behavior of "ScalaConfig"

  val exampleConfig = ConfigFactory.parseString(
    """
      |{
      |  configVal: {}
      |  stringVal: "string"
      |  booleanVal: true
      |  intVal: 123
      |  longVal: 123
      |  doubleVal: 12.3
      |  bytesVal: 1G
      |}
    """.stripMargin)

  it should "return the config value as options when present" in {
    exampleConfig.getConfigOption("configVal") should be('defined)
    exampleConfig.getStringOption("stringVal") should be(Some("string"))
    exampleConfig.getBooleanOption("booleanVal") should be(Some(true))
    exampleConfig.getIntOption("intVal") should be(Some(123))
    exampleConfig.getLongOption("longVal") should be(Some(123L))
    exampleConfig.getDoubleOption("doubleVal") should be(Some(12.3))
    exampleConfig.getBytesOption("bytesVal") should be(Some(1024 * 1024 * 1024L))
  }

  it should "return the config as None when missing" in {
    exampleConfig.getConfigOption("missing") should be(None)
    exampleConfig.getStringOption("missing") should be(None)
    exampleConfig.getBooleanOption("missing") should be(None)
    exampleConfig.getIntOption("missing") should be(None)
    exampleConfig.getLongOption("missing") should be(None)
    exampleConfig.getDoubleOption("missing") should be(None)
    exampleConfig.getBytesOption("missing") should be(None)
  }

  it should "return the config value when found" in {
    exampleConfig.getConfigOr("configVal") should be(ScalaConfig.empty)
    exampleConfig.getStringOr("stringVal") should be("string")
    exampleConfig.getBooleanOr("booleanVal") should be(right = true)
    exampleConfig.getIntOr("intVal") should be(123)
    exampleConfig.getLongOr("longVal") should be(123L)
    exampleConfig.getDoubleOr("doubleVal") should be(12.3)
    exampleConfig.getBytesOr("bytesVal") should be(1024 * 1024 * 1024L)
  }

  it should "not execute the default when the config value is found" in {
    exampleConfig.getConfigOr("configVal", fail()) should be(ScalaConfig.empty)
    exampleConfig.getStringOr("stringVal", fail()) should be("string")
    exampleConfig.getBooleanOr("booleanVal", fail()) should be(right = true)
    exampleConfig.getIntOr("intVal", fail()) should be(123)
    exampleConfig.getLongOr("longVal", fail()) should be(123L)
    exampleConfig.getDoubleOr("doubleVal", fail()) should be(12.3)
    exampleConfig.getBytesOr("bytesVal", fail()) should be(1024 * 1024 * 1024L)
  }

  it should "return the implicit defaults when missing" in {
    exampleConfig.getConfigOr("missing") should be(ScalaConfig.empty)
    exampleConfig.getStringOr("missing") should be("")
    exampleConfig.getBooleanOr("missing") should be(right = false)
    exampleConfig.getIntOr("missing") should be(0)
    exampleConfig.getLongOr("missing") should be(0L)
    exampleConfig.getDoubleOr("missing") should be(0.0)
    exampleConfig.getBytesOr("missing") should be(0L)
  }

  it should "return the explict defaults when missing" in {
    exampleConfig.getConfigOr("missing", ConfigFactory.parseString("""hello: world""")).
      getString("hello") should be("world")
    exampleConfig.getStringOr("missing", "default") should be("default")
    exampleConfig.getBooleanOr("missing", default = true) should be(right = true)
    exampleConfig.getIntOr("missing", 321) should be(321)
    exampleConfig.getLongOr("missing", 321L) should be(321L)
    exampleConfig.getDoubleOr("missing", 32.1) should be(32.1)
    exampleConfig.getBytesOr("missing", 321L) should be(321L)
  }

  it should "be ready to throw an unexpected getOption error" in {
    class ExpectedException extends RuntimeException
    intercept[ExpectedException](ScalaConfig.getOption("not used", (_: String) => throw new ExpectedException))
  }
}
