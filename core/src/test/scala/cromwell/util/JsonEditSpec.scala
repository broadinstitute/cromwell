package cromwell.util

import JsonEditor._
import cats.data.NonEmptyList
import io.circe.{Json, ParsingFailure}
import io.circe.parser._
import org.scalatest.{FlatSpec, Matchers}
import cats.syntax.either._

class JsonEditSpec extends FlatSpec with Matchers{
  import JsonEditSpec._

  def testJson(f: Json => Json): Either[String, Json] =
    for {
      json <- contrivedJsonEither
      newJson = f(json)
    }  yield newJson

  def testJsonAndGetKeys(f: Json => Json): Either[String, Iterable[String]] = {
    for {
      newJson <- testJson(f)
      keys <- newJson.hcursor.keys.toRight("no keys found!")
    }  yield keys
  }


  "Json Munger" should "remove excludes" in {
      val either = testJsonAndGetKeys(includeExcludeJson(_, None, Some(NonEmptyList.one("foo"))))
      assert(either.right.get.head === "other")
    }

  it should "remove nested keys excludes" in {
    val either = testJson(excludeJson(_, NonEmptyList.one("deep")))
    assert(either.right.get.hcursor.downField("nested").downField("inner").keys.get.head === "keepme")
  }

  it should "remove multiple nested keys excludes" in {
    val either = testJson(excludeJson(_, NonEmptyList.of("deep", "wildcard")))
    val keys = either.right.get.hcursor.downField("nested").downField("inner").keys.get
    assert(keys.head === "keepme")
    assert(keys.size === 1)
  }

  it should "keep includes" in {
    val either = testJsonAndGetKeys(includeExcludeJson(_, Some(NonEmptyList.one("foo")), None))
    assert(either.right.get.head === "foo")
  }

  it should "keep nested includes" in {
    val either = testJson(includeJson(_, NonEmptyList.one("keepme")))
    assert(either.right.get.hcursor.downField("nested").downField("inner").keys.get.head === "keepme")
  }

  it should "keep multiple nested includes" in {
    val either = testJson(includeJson(_, NonEmptyList.of("keepme", "wildcard")))
    val keys = either.right.get.hcursor.downField("nested").downField("inner").keys.get
    assert(keys.head === "keepme")
    assert(keys.tail.head === "wildcard")
  }

  it should "keep outputs only" in {
    val cursor = parse(helloWorldJsonOutput).map(outputs).right.get.hcursor
    val keys = cursor.keys
    println(s"keys $keys")
    assert(keys.get.size === 2)
    assert(keys.get.toSet.contains("outputs") === true)
    assert(keys.get.toSet.contains("id") === true)
    assert(cursor.downField("outputs").keys.get.head === "main_workflow.main_output")
  }

  it should "add labels" in {
    val newJson = realJson.map(augmentLabels(_, Map(("new","label")))).right.get
    val newLabels = newJson.hcursor.downField("labels").keys.get
    assert(newLabels.size  === 2)
  }
}

object JsonEditSpec {

  val contrivedJson =
    """
      {
        "foo": "bar",
         "other":"baz",
         "nested": {
           "inner": {
             "deep":"some",
             "keepme": "more",
             "wildcard": "again"
           }
         }
       }
      """.stripMargin

  val contrivedJsonEither: Either[String, Json] = parse(contrivedJson).leftMap(_.toString)

  val helloWorldJsonOutput =
    """
      |{
      |  "workflowName": "main_workflow",
      |  "actualWorkflowLanguageVersion": "draft-2",
      |  "submittedFiles": {
      |    "workflow": "import \"sub_workflow_hello_world_import.wdl\" as sub\n\nworkflow main_workflow {\n    call sub.wf_hello { input: wf_hello_input = \"sub world\" }\n    output {\n        String main_output = wf_hello.salutation\n    }\n}",
      |    "root": "",
      |    "options": "{\n\n}",
      |    "inputs": "{}",
      |    "workflowUrl": "",
      |    "labels": "{}",
      |    "imports": {
      |      "sub_workflow_hello_world_import.wdl": "task hello {\n  String addressee\n  command {\n    echo \"Hello ${addressee}!\"\n  }\n  runtime {\n      docker: \"ubuntu:latest\"\n  }\n  output {\n    String salutation = read_string(stdout())\n  }\n}\n\nworkflow wf_hello {\n  String wf_hello_input = \"world\"\n  \n  call hello { input: addressee = wf_hello_input }\n  \n  output {\n    String salutation = hello.salutation\n  }\n}\n"
      |    }
      |  },
      |  "calls": {
      |    "main_workflow.wf_hello": [
      |      {
      |        "executionStatus": "Done",
      |        "shardIndex": -1,
      |        "outputs": {
      |          "salutation": "Hello sub world!"
      |        },
      |        "inputs": {
      |          "wf_hello_input": "sub world"
      |        },
      |        "end": "2019-07-22T13:32:39.610-04:00",
      |        "attempt": 1,
      |        "executionEvents": [
      |          {
      |            "startTime": "2019-07-22T13:32:24.140-04:00",
      |            "endTime": "2019-07-22T13:32:24.147-04:00",
      |            "description": "WaitingForValueStore"
      |          },
      |          {
      |            "startTime": "2019-07-22T13:32:24.147-04:00",
      |            "endTime": "2019-07-22T13:32:24.196-04:00",
      |            "description": "SubWorkflowPreparingState"
      |          },
      |          {
      |            "startTime": "2019-07-22T13:32:24.134-04:00",
      |            "description": "SubWorkflowPendingState",
      |            "endTime": "2019-07-22T13:32:24.140-04:00"
      |          },
      |          {
      |            "startTime": "2019-07-22T13:32:24.196-04:00",
      |            "description": "SubWorkflowRunningState",
      |            "endTime": "2019-07-22T13:32:39.609-04:00"
      |          }
      |        ],
      |        "start": "2019-07-22T13:32:24.133-04:00",
      |        "subWorkflowId": "ba4a8cec-889e-4cce-9f89-d7df2a12d8d4"
      |      }
      |    ]
      |  },
      |  "outputs": {
      |    "main_workflow.main_output": "Hello sub world!"
      |  },
      |  "workflowRoot": "/home/dan/cromwell/cromwell-executions/main_workflow/757d0bcc-b636-4658-99d4-9b4b3767f1d1",
      |  "actualWorkflowLanguage": "WDL",
      |  "id": "757d0bcc-b636-4658-99d4-9b4b3767f1d1",
      |  "inputs": {},
      |  "labels": {
      |    "cromwell-workflow-id": "cromwell-757d0bcc-b636-4658-99d4-9b4b3767f1d1"
      |  },
      |  "submission": "2019-07-22T13:32:02.123-04:00",
      |  "status": "Succeeded",
      |  "end": "2019-07-22T13:32:41.529-04:00",
      |  "start": "2019-07-22T13:32:20.434-04:00"
      |}
    """.stripMargin

  val realJson: Either[ParsingFailure, Json] = parse(helloWorldJsonOutput)
}
