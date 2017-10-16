package wdl.util

import org.scalatest.{FlatSpec, Matchers}
import wom.types._
import wom.util.JsUtil
import wom.values._

class JsUtilSpec extends FlatSpec with Matchers {

  behavior of "JsUtil"

  it should "eval" in {
    val values = Map(
      "myName" -> WomMap(
        WomMapType(WomBooleanType, WomArrayType(WomStringType)),
        Map(WomBoolean(true) -> WomArray(WomArrayType(WomStringType), Seq(WomString("myValue"))))
      )
    )

    val expr = "myName[true][0] + 'Plus'"

    val result: WomValue = JsUtil.eval(expr, values)

    result should be(WomString("myValuePlus"))
  }

}
