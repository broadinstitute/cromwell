package cwl

import org.scalatest.{FlatSpec, Matchers}
import wom.values._

import collection.JavaConverters._

class JsUtilSpec extends FlatSpec with Matchers {

  behavior of "JsUtil"

  it should "eval" in {
    val values = Map(
      "myName" -> Map(true ->  Array("myValue")).asJava.asInstanceOf[AnyRef]
    ).asJava

    val expr = "myName[true][0] + 'Plus'"

    val result: WomValue = JsUtil.eval(expr, values)

    result should be(WomString("myValuePlus"))
  }

  it should "eval inputs of different types" in {
    val values: Map[String, WomValue] = Map(
      "myName" -> WomString("hi"),
      "someother" -> WomBoolean(false)
    )

    ParameterContext().addInputs(values).foreach{pc =>

      val expr = "inputs.myName"

      val result: WomValue = JsUtil.eval(expr, pc.ecmaScriptValues)

      result should be(WomString("hi"))

      val expr2 = "inputs.someother"

      val result2: WomValue = JsUtil.eval(expr2, pc.ecmaScriptValues)

      result2 should be(WomBoolean(false))
    }

  }

}
