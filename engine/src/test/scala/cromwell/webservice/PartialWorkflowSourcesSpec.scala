package cromwell.webservice

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import spray.json.DefaultJsonProtocol._
import spray.json._

class PartialWorkflowSourcesSpec extends FlatSpec with Matchers {

  behavior of "mergeMaps method"
  it should "successfully merge and override multiple input files" in {
    val input1 = Map("wf.a1" -> "hello", "wf.a2" -> "world").toJson.toString
    val input2 = Map.empty[String, String].toJson.toString
    val overrideInput1 = Map("wf.a2" -> "universe").toJson.toString
    val mergedMapsErrorOr = PartialWorkflowSources.mergeMaps(Seq(Option(input1), Option(input2), Option(overrideInput1)))

    mergedMapsErrorOr match {
      case Valid(inputs) => {
        inputs.fields.keys should contain allOf("wf.a1", "wf.a2")
        inputs.fields("wf.a2") should be(JsString("universe"))
      }
      case Invalid(error) => fail(s"This is unexpected! This test should pass! Error: $error")
    }
  }

  it should "return error when workflow input is not a valid json object" in {
    val invalidJsonInput = "\"invalidInput\""
    val mergedMapsErrorOr = PartialWorkflowSources.mergeMaps(Seq(Option(invalidJsonInput)))

    mergedMapsErrorOr match {
      case Valid(_) => fail("This is unexpected! This test is designed to fail!")
      case Invalid(error) => error.head shouldBe "Submitted input '\"invalidInput\"' of type JsString is not a JSON object."
    }
  }

  it should "return error when auxiliary workflow input is not a valid json" in {
    val validInput = "{\"key\":\"value\"}"
    val invalidAuxJsonInput = "invalidInput"
    val mergedMapsErrorOr = PartialWorkflowSources.mergeMaps(Seq(Option(validInput), Option(invalidAuxJsonInput)))

    mergedMapsErrorOr match {
      case Valid(_) => fail("This is unexpected! This test is designed to fail!")
      case Invalid(error) => error.head shouldBe "Failed to parse input: 'invalidInput', which is not a valid json. Please check for syntactical errors. (reason 1 of 1): Unexpected character 'i' at input index 0 (line 1, position 1), expected JSON Value:\ninvalidInput\n^\n"
    }
  }
}
