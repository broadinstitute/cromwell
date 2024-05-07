package cromwell.webservice

import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import common.assertion.CromwellTimeoutSpec
import cromwell.core.WorkflowId
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json.DefaultJsonProtocol._
import spray.json._

class PartialWorkflowSourcesSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "mergeMaps method"
  it should "successfully merge and override multiple input files" in {
    val input1 = Map("wf.a1" -> "hello", "wf.a2" -> "world").toJson.toString
    val input2 = Map.empty[String, String].toJson.toString
    val overrideInput1 = Map("wf.a2" -> "universe").toJson.toString
    val mergedMapsErrorOr =
      PartialWorkflowSources.mergeMaps(Seq(Option(input1), Option(input2), Option(overrideInput1)))

    mergedMapsErrorOr match {
      case Valid(inputs) =>
        inputs.fields.keys should contain allOf ("wf.a1", "wf.a2")
        inputs.fields("wf.a2") should be(JsString("universe"))
      case Invalid(error) => fail(s"This is unexpected! This test should pass! Error: $error")
    }
  }

  it should "return error when workflow input is not a valid json object" in {
    val invalidJsonInput = "\"invalidInput\""
    val mergedMapsErrorOr = PartialWorkflowSources.mergeMaps(Seq(Option(invalidJsonInput)))

    mergedMapsErrorOr match {
      case Valid(_) => fail("This is unexpected! This test is designed to fail!")
      case Invalid(error) =>
        error.head shouldBe "Submitted input '\"invalidInput\"' of type JsString is not a valid JSON object."
    }
  }

  it should "return error when auxiliary workflow input is not a valid json" in {
    val validInput = "{\"key\":\"value\"}"
    val invalidAuxJsonInput = "invalidInput"
    val mergedMapsErrorOr = PartialWorkflowSources.mergeMaps(Seq(Option(validInput), Option(invalidAuxJsonInput)))

    mergedMapsErrorOr match {
      case Valid(_) => fail("This is unexpected! This test is designed to fail!")
      case Invalid(error) =>
        error.head shouldBe "Failed to parse input: 'invalidInput', which is not a valid json. Please check for syntactical errors. (reason 1 of 1): Unexpected character 'i' at input index 0 (line 1, position 1), expected JSON Value:\ninvalidInput\n^\n"
    }
  }

  behavior of "workflow inputs parsing"

  it should "consider a single json object as a single workflow inputs object" in {
    val input =
      """{
        |  "mywf.inInt": 1,
        |  "mywf.inString": "one"
        |}
        |""".stripMargin

    val expected = Vector("""{"mywf.inInt":1,"mywf.inString":"one"}""").validNel

    val actual = PartialWorkflowSources.workflowInputsValidation(input)

    actual should be(expected)
  }

  it should "consider an array of json objects as a list of workflow inputs objects" in {
    val input =
      """[{
        |  "mywf.inInt": 1,
        |  "mywf.inString": "one"
        |}, {
        |  "mywf.inInt": 2,
        |  "mywf.inString": "two"
        |}]
        |""".stripMargin

    val expected =
      Vector("""{"mywf.inInt":1,"mywf.inString":"one"}""", """{"mywf.inInt":2,"mywf.inString":"two"}""").validNel

    val actual = PartialWorkflowSources.workflowInputsValidation(input)

    actual should be(expected)
  }

  behavior of "requested workflow ID parsing"

  it should "interpret a workflow ID string" in {
    val input = """ "31db02f4-82b8-11ec-a8a3-0242ac120002" """

    val expected = Vector(WorkflowId.fromString("31db02f4-82b8-11ec-a8a3-0242ac120002")).validNel

    val actual = PartialWorkflowSources.workflowIdsValidation(input)

    actual should be(expected)
  }

  it should "interpret a workflow ID array" in {
    val input =
      """ [
        |"57adac0c-6c02-4035-b91e-75fe217a5662",
        |"5fd418c6-42ac-47a3-bb59-d97fd0a3e864",
        |"ca808e1e-3227-4315-b713-cbfa8e598e74"
        | ]""".stripMargin

    val expected = Vector(
      WorkflowId.fromString("57adac0c-6c02-4035-b91e-75fe217a5662"),
      WorkflowId.fromString("5fd418c6-42ac-47a3-bb59-d97fd0a3e864"),
      WorkflowId.fromString("ca808e1e-3227-4315-b713-cbfa8e598e74")
    ).validNel

    val actual = PartialWorkflowSources.workflowIdsValidation(input)

    actual should be(expected)
  }
}
