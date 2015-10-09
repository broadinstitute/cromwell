package lenthall.util

import com.typesafe.config.ConfigFactory
import lenthall.util.ScalaConfig._
import org.scalatest.{FlatSpec, Matchers}

class ScalaConfigSpec extends FlatSpec with Matchers {
  behavior of "ScalaConfig"

  val exampleConfig = ConfigFactory.parseString(
    """
      |{
      |  stringVal: "string"
      |  booleanVal: true
      |  intVal: 123
      |  longVal: 123
      |  doubleVal: 12.3
      |}
    """.stripMargin)

  it should "return the config value as options when present" in {
    exampleConfig.getStringOption("stringVal") should be(Some("string"))
    exampleConfig.getBooleanOption("booleanVal") should be(Some(true))
    exampleConfig.getIntOption("intVal") should be(Some(123))
    exampleConfig.getLongOption("longVal") should be(Some(123L))
    exampleConfig.getDoubleOption("doubleVal") should be(Some(12.3))
  }

  it should "return the config as None when missing" in {
    exampleConfig.getStringOption("missing") should be(None)
    exampleConfig.getBooleanOption("missing") should be(None)
    exampleConfig.getIntOption("missing") should be(None)
    exampleConfig.getLongOption("missing") should be(None)
    exampleConfig.getDoubleOption("missing") should be(None)
  }

  it should "return the config value when found" in {
    exampleConfig.getStringOr("stringVal") should be("string")
    exampleConfig.getBooleanOr("booleanVal") should be(right = true)
    exampleConfig.getIntOr("intVal") should be(123)
    exampleConfig.getLongOr("longVal") should be(123L)
    exampleConfig.getDoubleOr("doubleVal") should be(12.3)
  }

  it should "return the implict defaults when missing" in {
    exampleConfig.getStringOr("missing") should be("")
    exampleConfig.getBooleanOr("missing") should be(right = false)
    exampleConfig.getIntOr("missing") should be(0)
    exampleConfig.getLongOr("missing") should be(0L)
    exampleConfig.getDoubleOr("missing") should be(0.0)
  }

  it should "return the explict defaults when missing" in {
    exampleConfig.getStringOr("missing", "default") should be("default")
    exampleConfig.getBooleanOr("missing", default = true) should be(right = true)
    exampleConfig.getIntOr("missing", 321) should be(321)
    exampleConfig.getLongOr("missing", 321L) should be(321L)
    exampleConfig.getDoubleOr("missing", 32.1) should be(32.1)
  }
}
