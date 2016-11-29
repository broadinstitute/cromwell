package cromwell.util

import scala.Vector

import org.scalatest.FlatSpec
import org.scalatest.Matchers

import JsonFormatting.WdlValueJsonFormatter.WdlValueJsonFormat
import spray.json.{ JsObject, pimpString }
import wdl4s.types.{ WdlArrayType, WdlStringType }
import wdl4s.values.{ WdlArray, WdlPair, WdlString }

class WdlValueJsonFormatterSpec extends FlatSpec with Matchers {

  behavior of "WdlValueJsonFormat"

  it should "write WdlPair to left/right structured JsObject" in {
    val left = "sanders"
    val right = Vector("rubio", "carson", "cruz")
    val wdlPair = WdlPair(WdlString(left), WdlArray(WdlArrayType(WdlStringType), right.map { WdlString(_) }))
    val ExpectedJson: JsObject =
      """|{
         |  "left": "sanders",
         |  "right": ["rubio", "carson", "cruz"]
         |}""".stripMargin.parseJson.asJsObject
    WdlValueJsonFormat.write(wdlPair) should matchPattern { case ExpectedJson => }
  }
}
