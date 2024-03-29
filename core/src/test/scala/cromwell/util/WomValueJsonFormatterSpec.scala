package cromwell.util

import common.assertion.CromwellTimeoutSpec
import cromwell.util.JsonFormatting.WomValueJsonFormatter.WomValueJsonFormat
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json.{enrichString, JsObject}
import wom.types._
import wom.values._

class WomValueJsonFormatterSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "WdlValueJsonFormat"

  it should "write WdlPair to left/right structured JsObject" in {
    val left = "sanders"
    val right = Vector("rubio", "carson", "cruz")
    val wdlPair = WomPair(WomString(left), WomArray(WomArrayType(WomStringType), right.map(WomString(_))))
    val ExpectedJson: JsObject =
      """|{
         |  "left": "sanders",
         |  "right": ["rubio", "carson", "cruz"]
         |}""".stripMargin.parseJson.asJsObject
    WomValueJsonFormat.write(wdlPair) should matchPattern { case ExpectedJson => }
  }
}
