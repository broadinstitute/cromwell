package cwl.internal

import common.validation.Validation._
import org.scalatest.{FlatSpec, Matchers}
import wom.types._
import wom.values._

class JsUtilSpec extends FlatSpec with Matchers {

  behavior of "JsUtil"

  it should "eval" in {
    val values =
      "myName" -> WomMap(
        WomMapType(WomBooleanType, WomArrayType(WomStringType)),
        Map(WomBoolean(true) -> WomArray(WomArrayType(WomStringType), Seq(WomString("myValue"))))
      )

    // NOTE: By using JsMap a JSObject, myName[true] doesn't work.
    val expr = """myName["true"][0] + 'Plus'"""

    val result: WomValue = JsUtil.evalStructish(expr,  values).toTry.get

    result should be(WomString("myValuePlus"))
  }

   it should "JSON.stringify" in {
    val values =
      "myName" -> WomMap(
        WomMapType(WomBooleanType, WomArrayType(WomStringType)),
        Map(WomBoolean(true) -> WomArray(WomArrayType(WomStringType), Seq(WomString("myValue"))))
      )



    val expr = "JSON.stringify(myName)"

    val result: WomValue = JsUtil.evalStructish(expr, values).toTry.get

    result should be(WomString("""{"true":["myValue"]}"""))
  }
}
