package wom.util

import common.validation.Validation._
import org.scalatest.{FlatSpec, Matchers}
import wom.types._
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

    // NOTE: By using JsMap a JSObject, myName[true] doesn't work.
    val expr = """myName["true"][0] + 'Plus'"""

    val result: WomValue = JsUtil.eval(expr, values).toTry.get

    result should be(WomString("myValuePlus"))
  }

  it should "eval arrays" in {
    val values = Map(
      "myName" -> WomArray(Seq(WomInteger(1),WomInteger(3),WomInteger(2)) )
    )

    val expr = """myName.sorted"""

    val result: WomValue = JsUtil.eval(expr, values).toTry.get

    result should be(WomArray(Seq(WomInteger(1),WomInteger(2),WomInteger(3))))
  }
}
