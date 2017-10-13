package wdl.util

import org.scalatest.{FlatSpec, Matchers}
import wom.types._
import wom.util.JsUtil
import wom.values._

class JsUtilSpec extends FlatSpec with Matchers {

  behavior of "JsUtil"

  it should "eval" in {
    val values = Map(
      "myName" -> WdlMap(
        WdlMapType(WdlBooleanType, WdlArrayType(WdlStringType)),
        Map(WdlBoolean(true) -> WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("myValue"))))
      )
    )

    val expr = "myName[true][0] + 'Plus'"

    val result: WdlValue = JsUtil.eval(expr, values)

    result should be(WdlString("myValuePlus"))
  }

}
