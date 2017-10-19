package cromwell.util

import cromwell.util.JsonFormatting.WdlValueJsonFormatter.WdlValueJsonFormat
import org.scalatest.{FlatSpec, Matchers}
import spray.json.{JsObject, pimpString}
import wom.types._
import wom.values._

class WomValueJsonFormatterSpec extends FlatSpec with Matchers {

  behavior of "WdlValueJsonFormat"

  it should "write WdlPair to left/right structured JsObject" in {
    val left = "sanders"
    val right = Vector("rubio", "carson", "cruz")
    val wdlPair = WomPair(WomString(left), WomArray(WomArrayType(WomStringType), right.map { WomString(_) }))
    val ExpectedJson: JsObject =
      """|{
         |  "left": "sanders",
         |  "right": ["rubio", "carson", "cruz"]
         |}""".stripMargin.parseJson.asJsObject
    WdlValueJsonFormat.write(wdlPair) should matchPattern { case ExpectedJson => }
  }
}
