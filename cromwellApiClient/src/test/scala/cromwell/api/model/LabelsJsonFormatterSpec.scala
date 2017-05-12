package cromwell.api

import cromwell.api.model._
import spray.json._
import org.scalatest.FlatSpec
import org.scalatest.Matchers

class LabelsJsonFormatterSpec extends FlatSpec with Matchers {
  import cromwell.api.model.LabelsJsonFormatter._

  behavior of "WdlValueJsonFormat"

  it should "write a Label as a structured JsObject" in {
    val label = List(Label("test-key", "test-value"))
    val expectedJson: JsObject =
      """|{
         |  "test-key": "test-value"
         |}""".stripMargin.parseJson.asJsObject

    label.toJson shouldEqual expectedJson
  }

  it should "write an optional list of labels as a structured JsObject" in {
    val labels =  Option(List(Label("key-1", "value-1"),
                              Label("key-2", "value-2"),
                              Label("key-3", "value-3")))

    val expectedJson: JsObject =
      """|{
        |  "key-1":"value-1",
        |  "key-2":"value-2",
        |  "key-3":"value-3"
        |}""".stripMargin.parseJson.asJsObject

    labels.toJson shouldEqual expectedJson
  }

  it should "read a list of labels as a structured JsObject" in {
    val expectedLabels =  List(Label("key-1", "value-1"),
      Label("key-2", "value-2"),
      Label("key-3", "value-3"))

    val labelsJson: JsObject =
      """|{
        |  "key-1":"value-1",
        |  "key-2":"value-2",
        |  "key-3":"value-3"
        |}""".stripMargin.parseJson.asJsObject

    labelsJson.convertTo[List[Label]] shouldBe expectedLabels
  }
}