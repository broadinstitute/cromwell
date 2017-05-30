package cromwell.api.model

import spray.json._
import org.scalatest.FlatSpec
import org.scalatest.Matchers

class LabelsJsonFormatterSpec extends FlatSpec with Matchers {
  import cromwell.api.model.LabelsJsonFormatter._

  behavior of "WdlValueJsonFormat"

  val sampleLabels = List(Label("key-1", "value-1"), Label("key-2", "value-2"), Label("key-3", "value-3"))
  val sampleJson =  """|{
                       |  "key-1":"value-1",
                       |  "key-2":"value-2",
                       |  "key-3":"value-3"
                       |}""".stripMargin.parseJson.asJsObject

  it should "write a Label as a structured JsObject" in {
    val label = List(Label("test-key", "test-value"))
    val expectedJson: JsObject =
      """|{
         |  "test-key": "test-value"
         |}""".stripMargin.parseJson.asJsObject

    label.toJson shouldEqual expectedJson
  }

  it should "write an optional list of labels as a structured JsObject" in {
    Option(sampleLabels).toJson shouldEqual sampleJson
  }

  it should "read a list of labels as a structured JsObject" in {
    sampleJson.convertTo[List[Label]] shouldBe sampleLabels
  }
}
