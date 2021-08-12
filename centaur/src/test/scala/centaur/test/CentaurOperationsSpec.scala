package centaur.test

import java.util.UUID

import centaur.test.workflow.Workflow
import cromwell.api.model.{SubmittedWorkflow, WorkflowId}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json._

import scala.util.{Failure, Success}
import scala.concurrent.Await
import scala.concurrent.duration._

class CentaurOperationsSpec extends AnyFlatSpec with Matchers {
  behavior of "validateMetadataJson"

  val placeholderSubmittedWorkflow: SubmittedWorkflow = SubmittedWorkflow(id = WorkflowId(UUID.randomUUID()), null, null)
  val placeholderWorkflow: Workflow = Workflow(testName = "", null, null, null, null, null, false, false, false, null)

  val allowableOneWordAdditions = List("farmer")

  def runTest(json1: String, json2: String, expectMatching: Boolean): Unit = {
    val validation = Operations.validateMetadataJson("",
      json1.parseJson.asJsObject,
      json2.parseJson.asJsObject,
      placeholderSubmittedWorkflow,
      placeholderWorkflow,
      allowableAddedOneWordFields = allowableOneWordAdditions).unsafeToFuture()
    Await.ready(validation, atMost = 10.seconds)
    validation.value.get match {
      case Success(()) if expectMatching => // great
      case Success(_) if !expectMatching => fail("Metadata unexpectedly matches")
      case Failure(e) if expectMatching  => fail("Metadata unexpectedly mismatches", e)
      case Failure(_) if !expectMatching => // great
    }
  }

  it should "decide identical JSON as matching" in {
    val json1 =
      """{
        |  "id": "foobar",
        |  "calls": {
        |    "workflow.callname": [
        |      {
        |        "callfield1": "x",
        |        "callfield2": "y"
        |      }
        |    ]
        |  }
        |}""".stripMargin

    val json2 = json1
    runTest(json1, json2, expectMatching = true)
  }

  it should "decide different JSON as non-matching" in {
    val json1 =
      """{
        |  "id": "foobar",
        |  "calls": {
        |    "workflow.callname": [
        |      {
        |        "callfield1": "x",
        |        "callfield2": "y"
        |      }
        |    ]
        |  }
        |}""".stripMargin

    val json2 =
      """{
        |  "id": "foobar",
        |  "calls": {
        |    "workflow.callname": [
        |      {
        |        "callfield1": "x",
        |        "callfield2": "y",
        |        "foo": "this is bad"
        |      }
        |    ]
        |  },
        |  "foo": "this is also bad"
        |}""".stripMargin

    runTest(json1, json2, expectMatching = false)
  }

  it should "ignore allowable additions at the root level and within calls" in {
    val json1 =
      """{
        |  "id": "foobar",
        |  "calls": {
        |    "workflow.callname": [
        |      {
        |        "callfield1": "x",
        |        "callfield2": "y"
        |      }
        |    ]
        |  }
        |}""".stripMargin

    val json2 =
      """{
        |  "id": "foobar",
        |  "calls": {
        |    "workflow.callname": [
        |      {
        |        "callfield1": "x",
        |        "callfield2": "y",
        |        "farmer": "hoggett"
        |      }
        |    ]
        |  },
        |  "farmer": "hoggett"
        |}""".stripMargin

    runTest(json1, json2, expectMatching = true)
  }
}
