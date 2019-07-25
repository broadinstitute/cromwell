package cromwell.util

import JsonEditor._
import cats.data.NonEmptyList
import io.circe.{Json, ParsingFailure}
import io.circe.parser._
import org.scalatest.{FlatSpec, Matchers}
import cats.syntax.either._

class JsonEditorSpec extends FlatSpec with Matchers{
  import JsonEditorSpec._

  def testJson(f: Json => Json): Either[String, Json] =
    for {
      json <- contrivedJsonEither
      newJson = f(json)
    }  yield newJson

  def testJsonAndGetKeys(f: Json => Json): Either[String, Iterable[String]] =
    for {
      newJson <- testJson(f)
      keys <- newJson.hcursor.keys.toRight("no keys found!")
    }  yield keys

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

  it should "remove subworkflow info" in {
    println(subWorkflowJson.map(removeSubworkflowData))
    val sub = subWorkflowJson.map(removeSubworkflowData).right.get
    val keys = sub.hcursor.downField("calls").downField("sub_workflow_interactions.countEvens").downArray.keys
    assert(keys.contains("subWorkflowMetadata") === false)
  }
}

object JsonEditorSpec {

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

  val metadataWithSubWorkflowMetadata = """{
                                          |  "workflowName": "sub_workflow_interactions",
                                          |  "workflowProcessingEvents": [
                                          |    {
                                          |      "cromwellId": "cromid-c58ca80",
                                          |      "description": "Finished",
                                          |      "cromwellVersion": "45-9409587-SNAP",
                                          |      "timestamp": "2019-07-23T07:43:43.712Z"
                                          |    },
                                          |    {
                                          |      "cromwellId": "cromid-c58ca80",
                                          |      "description": "PickedUp",
                                          |      "timestamp": "2019-07-23T07:42:46.342Z",
                                          |      "cromwellVersion": "45-9409587-SNAP"
                                          |    }
                                          |  ],
                                          |  "actualWorkflowLanguageVersion": "draft-2",
                                          |  "submittedFiles": {
                                          |    "workflow": "import \"sub_workflow_interactions_import.wdl\" as counter\n\ntask hello {\n  String addressee\n  \n  command {\n    echo \"Hello ${addressee}!\" > hello\n    wc -w < hello > count\n  }\n  runtime {\n      docker: \"ubuntu:latest\"\n  }\n  output {\n    String salutation = read_string(\"hello\")\n    Int count = read_int(\"count\")\n  }\n}\n\nworkflow sub_workflow_interactions {\n  call hello { input: addressee = \"Sub Workflow World\" }\n  call counter.countEvens { input: max = hello.count } # Sub workflow depends on previous task call\n  call hello as secondHello { input: addressee = countEvens.someStringOutput } # Task call depends on previous sub workflow call\n  \n  output {\n    # old output syntax\n    hello.* \n    \n    # new output syntax\n    String out2 = secondHello.salutation\n    Array[Int] out4 = read_lines(countEvens.evenFile)\n  }\n}\n",
                                          |    "root": "",
                                          |    "options": "{\n\n}",
                                          |    "inputs": "{}",
                                          |    "workflowUrl": "https://raw.githubusercontent.com/broadinstitute/cromwell/develop/centaur/src/main/resources/standardTestCases/sub_workflow_interactions/sub_workflow_interactions.wdl",
                                          |    "labels": "{}",
                                          |    "imports": {
                                          |      "sub_workflow_interactions_import.wdl": "task countTo {\n    Int value\n    command {\n        seq 0 1 ${value}\n    }\n    runtime {\n          docker: \"ubuntu:latest\"\n      }\n    output {\n        File range = stdout()\n    }\n}\n\ntask filterEvens {\n    File numbers\n    command {\n        grep '[02468]$' ${numbers} > evens\n    }\n    runtime {\n          docker: \"ubuntu:latest\"\n      }\n    output {\n        File evens = \"evens\"\n    }\n}\n\nworkflow countEvens {\n    Int max = 10\n    \n    call countTo { input: value = max }\n    call filterEvens { input: numbers = countTo.range }\n    output {\n        String someStringOutput = \"I'm an output\"\n        File evenFile = filterEvens.evens\n    }\n}\n"
                                          |    }
                                          |  },
                                          |  "calls": {
                                          |    "sub_workflow_interactions.secondHello": [
                                          |      {
                                          |        "executionStatus": "Done",
                                          |        "stdout": "/Users/gentrj/projects/cromwell/cromwell-executions/sub_workflow_interactions/2f070026-0186-4c8b-a577-59d818b5ac7c/call-secondHello/execution/stdout",
                                          |        "backendStatus": "Done",
                                          |        "compressedDockerSize": 26720871,
                                          |        "commandLine": "echo \"Hello I'm an output!\" > hello\nwc -w < hello > count",
                                          |        "shardIndex": -1,
                                          |        "outputs": {
                                          |          "count": 4,
                                          |          "salutation": "Hello I'm an output!"
                                          |        },
                                          |        "runtimeAttributes": {
                                          |          "maxRetries": "0",
                                          |          "continueOnReturnCode": "0",
                                          |          "docker": "ubuntu:latest",
                                          |          "failOnStderr": "false"
                                          |        },
                                          |        "callCaching": {
                                          |          "allowResultReuse": false,
                                          |          "effectiveCallCachingMode": "CallCachingOff"
                                          |        },
                                          |        "inputs": {
                                          |          "addressee": "I'm an output"
                                          |        },
                                          |        "returnCode": 0,
                                          |        "jobId": "17755",
                                          |        "backend": "Local",
                                          |        "end": "2019-07-23T07:43:42.326Z",
                                          |        "dockerImageUsed": "ubuntu@sha256:9b1702dcfe32c873a770a32cfd306dd7fc1c4fd134adfb783db68defc8894b3c",
                                          |        "stderr": "/Users/gentrj/projects/cromwell/cromwell-executions/sub_workflow_interactions/2f070026-0186-4c8b-a577-59d818b5ac7c/call-secondHello/execution/stderr",
                                          |        "callRoot": "/Users/gentrj/projects/cromwell/cromwell-executions/sub_workflow_interactions/2f070026-0186-4c8b-a577-59d818b5ac7c/call-secondHello",
                                          |        "attempt": 1,
                                          |        "executionEvents": [
                                          |          {
                                          |            "startTime": "2019-07-23T07:43:35.507Z",
                                          |            "description": "Pending",
                                          |            "endTime": "2019-07-23T07:43:35.508Z"
                                          |          },
                                          |          {
                                          |            "startTime": "2019-07-23T07:43:41.344Z",
                                          |            "description": "UpdatingJobStore",
                                          |            "endTime": "2019-07-23T07:43:42.326Z"
                                          |          },
                                          |          {
                                          |            "startTime": "2019-07-23T07:43:36.386Z",
                                          |            "description": "PreparingJob",
                                          |            "endTime": "2019-07-23T07:43:36.391Z"
                                          |          },
                                          |          {
                                          |            "startTime": "2019-07-23T07:43:36.386Z",
                                          |            "description": "WaitingForValueStore",
                                          |            "endTime": "2019-07-23T07:43:36.386Z"
                                          |          },
                                          |          {
                                          |            "startTime": "2019-07-23T07:43:35.508Z",
                                          |            "description": "RequestingExecutionToken",
                                          |            "endTime": "2019-07-23T07:43:36.386Z"
                                          |          },
                                          |          {
                                          |            "startTime": "2019-07-23T07:43:36.391Z",
                                          |            "description": "RunningJob",
                                          |            "endTime": "2019-07-23T07:43:41.344Z"
                                          |          }
                                          |        ],
                                          |        "start": "2019-07-23T07:43:35.507Z"
                                          |      }
                                          |    ],
                                          |    "sub_workflow_interactions.hello": [
                                          |      {
                                          |        "executionStatus": "Done",
                                          |        "stdout": "/Users/gentrj/projects/cromwell/cromwell-executions/sub_workflow_interactions/2f070026-0186-4c8b-a577-59d818b5ac7c/call-hello/execution/stdout",
                                          |        "backendStatus": "Done",
                                          |        "compressedDockerSize": 26720871,
                                          |        "commandLine": "echo \"Hello Sub Workflow World!\" > hello\nwc -w < hello > count",
                                          |        "shardIndex": -1,
                                          |        "outputs": {
                                          |          "salutation": "Hello Sub Workflow World!",
                                          |          "count": 4
                                          |        },
                                          |        "runtimeAttributes": {
                                          |          "docker": "ubuntu:latest",
                                          |          "maxRetries": "0",
                                          |          "failOnStderr": "false",
                                          |          "continueOnReturnCode": "0"
                                          |        },
                                          |        "callCaching": {
                                          |          "allowResultReuse": false,
                                          |          "effectiveCallCachingMode": "CallCachingOff"
                                          |        },
                                          |        "inputs": {
                                          |          "addressee": "Sub Workflow World"
                                          |        },
                                          |        "returnCode": 0,
                                          |        "jobId": "17716",
                                          |        "backend": "Local",
                                          |        "end": "2019-07-23T07:43:17.371Z",
                                          |        "dockerImageUsed": "ubuntu@sha256:9b1702dcfe32c873a770a32cfd306dd7fc1c4fd134adfb783db68defc8894b3c",
                                          |        "stderr": "/Users/gentrj/projects/cromwell/cromwell-executions/sub_workflow_interactions/2f070026-0186-4c8b-a577-59d818b5ac7c/call-hello/execution/stderr",
                                          |        "callRoot": "/Users/gentrj/projects/cromwell/cromwell-executions/sub_workflow_interactions/2f070026-0186-4c8b-a577-59d818b5ac7c/call-hello",
                                          |        "attempt": 1,
                                          |        "executionEvents": [
                                          |          {
                                          |            "endTime": "2019-07-23T07:42:51.381Z",
                                          |            "description": "RequestingExecutionToken",
                                          |            "startTime": "2019-07-23T07:42:50.512Z"
                                          |          },
                                          |          {
                                          |            "startTime": "2019-07-23T07:43:16.685Z",
                                          |            "description": "UpdatingJobStore",
                                          |            "endTime": "2019-07-23T07:43:17.372Z"
                                          |          },
                                          |          {
                                          |            "startTime": "2019-07-23T07:42:51.401Z",
                                          |            "description": "RunningJob",
                                          |            "endTime": "2019-07-23T07:43:16.685Z"
                                          |          },
                                          |          {
                                          |            "startTime": "2019-07-23T07:42:50.511Z",
                                          |            "description": "Pending",
                                          |            "endTime": "2019-07-23T07:42:50.512Z"
                                          |          },
                                          |          {
                                          |            "startTime": "2019-07-23T07:42:51.384Z",
                                          |            "description": "PreparingJob",
                                          |            "endTime": "2019-07-23T07:42:51.401Z"
                                          |          },
                                          |          {
                                          |            "description": "WaitingForValueStore",
                                          |            "startTime": "2019-07-23T07:42:51.381Z",
                                          |            "endTime": "2019-07-23T07:42:51.384Z"
                                          |          }
                                          |        ],
                                          |        "start": "2019-07-23T07:42:50.511Z"
                                          |      }
                                          |    ],
                                          |    "sub_workflow_interactions.countEvens": [
                                          |      {
                                          |        "executionStatus": "Done",
                                          |        "shardIndex": -1,
                                          |        "outputs": {
                                          |          "evenFile": "/Users/gentrj/projects/cromwell/cromwell-executions/sub_workflow_interactions/2f070026-0186-4c8b-a577-59d818b5ac7c/call-countEvens/counter.countEvens/b37615ef-9612-4f8c-95d8-d8ebffb5c287/call-filterEvens/execution/evens",
                                          |          "someStringOutput": "I'm an output"
                                          |        },
                                          |        "inputs": {
                                          |          "max": 4
                                          |        },
                                          |        "end": "2019-07-23T07:43:33.550Z",
                                          |        "attempt": 1,
                                          |        "executionEvents": [
                                          |          {
                                          |            "startTime": "2019-07-23T07:43:19.122Z",
                                          |            "description": "SubWorkflowPendingState",
                                          |            "endTime": "2019-07-23T07:43:19.124Z"
                                          |          },
                                          |          {
                                          |            "startTime": "2019-07-23T07:43:19.124Z",
                                          |            "description": "WaitingForValueStore",
                                          |            "endTime": "2019-07-23T07:43:19.128Z"
                                          |          },
                                          |          {
                                          |            "startTime": "2019-07-23T07:43:19.128Z",
                                          |            "description": "SubWorkflowPreparingState",
                                          |            "endTime": "2019-07-23T07:43:19.149Z"
                                          |          },
                                          |          {
                                          |            "startTime": "2019-07-23T07:43:19.149Z",
                                          |            "description": "SubWorkflowRunningState",
                                          |            "endTime": "2019-07-23T07:43:33.546Z"
                                          |          }
                                          |        ],
                                          |        "start": "2019-07-23T07:43:19.119Z",
                                          |        "subWorkflowId": "b37615ef-9612-4f8c-95d8-d8ebffb5c287"
                                          |      }
                                          |    ]
                                          |  },
                                          |  "outputs": {
                                          |    "sub_workflow_interactions.hello.salutation": "Hello Sub Workflow World!",
                                          |    "sub_workflow_interactions.hello.count": 4,
                                          |    "sub_workflow_interactions.out2": "Hello I'm an output!",
                                          |    "sub_workflow_interactions.out4": [
                                          |      0,
                                          |      2,
                                          |      4
                                          |    ]
                                          |  },
                                          |  "workflowRoot": "/Users/gentrj/projects/cromwell/cromwell-executions/sub_workflow_interactions/2f070026-0186-4c8b-a577-59d818b5ac7c",
                                          |  "actualWorkflowLanguage": "WDL",
                                          |  "id": "2f070026-0186-4c8b-a577-59d818b5ac7c",
                                          |  "inputs": {},
                                          |  "labels": {
                                          |    "cromwell-workflow-id": "cromwell-2f070026-0186-4c8b-a577-59d818b5ac7c"
                                          |  },
                                          |  "submission": "2019-07-23T07:42:29.195Z",
                                          |  "status": "Succeeded",
                                          |  "end": "2019-07-23T07:43:43.711Z",
                                          |  "start": "2019-07-23T07:42:46.352Z"
                                          |}
                                          |""".stripMargin

  val subWorkflowJson: Either[ParsingFailure, Json] = parse(metadataWithSubWorkflowMetadata)

}
