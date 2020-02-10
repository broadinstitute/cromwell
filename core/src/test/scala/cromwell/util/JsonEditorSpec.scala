package cromwell.util

import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import common.validation.ErrorOr.ErrorOr
import cromwell.util.JsonEditor._
import io.circe.parser._
import io.circe.{DecodingFailure, FailedCursor, Json}
import org.scalatest.{FlatSpec, Matchers}

import scala.io.Source

class JsonEditorSpec extends FlatSpec with Matchers {
  import JsonEditorSpec._

  "JsonEditor" should "remove excludes in workflows" in {
    val actual = excludeJson(helloWorldJson, NonEmptyList.of("calls", "submittedFiles", "actualWorkflowLanguage")).get
    val expectedMetadata =
      """
        |{
        |  "workflowName": "main_workflow",
        |  "actualWorkflowLanguageVersion": "draft-2",
        |  "outputs": {
        |    "main_workflow.main_output": "Hello sub world!"
        |  },
        |  "workflowRoot": "/home/dan/cromwell/cromwell-executions/main_workflow/757d0bcc-b636-4658-99d4-9b4b3767f1d1",
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
        |""".stripMargin
    val expected = parseString(expectedMetadata)
    actual shouldEqual expected
  }

  it should "remove nested excludes in workflows" in {
    // This test intentionally applies multiple colon separated filters with the same outermost term.
    val actual = excludeJson(helloWorldJson, NonEmptyList.of("calls", "submittedFiles:workflowUrl", "submittedFiles:imports:sub_workflow_hello_world_import.wdl", "actualWorkflowLanguage")).get
    val expectedMetadata =
      """
        |{
        |  "workflowName": "main_workflow",
        |  "actualWorkflowLanguageVersion": "draft-2",
        |  "submittedFiles": {
        |    "workflow": "import \"sub_workflow_hello_world_import.wdl\" as sub\n\nworkflow main_workflow {\n    call sub.wf_hello { input: wf_hello_input = \"sub world\" }\n    output {\n        String main_output = wf_hello.salutation\n    }\n}",
        |    "root": "",
        |    "options": "{\n\n}",
        |    "inputs": "{}",
        |    "labels": "{}"
        |  },
        |  "outputs": {
        |    "main_workflow.main_output": "Hello sub world!"
        |  },
        |  "workflowRoot": "/home/dan/cromwell/cromwell-executions/main_workflow/757d0bcc-b636-4658-99d4-9b4b3767f1d1",
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
        |
        |""".stripMargin
    val expected = parseString(expectedMetadata)
    actual shouldEqual expected
  }

  it should "remove nested excludes in workflows, calls and subworkflows" in {
    val actual = excludeJson(treblyNestedSubworkflowJson, NonEmptyList.of("callCaching:hashes:runtime attribute:docker", "executionStatus")).get
    val expected = parseMetadata("nested_excludes_trebly_nested_subworkflows.json")

    actual shouldEqual expected
  }

  it should "remove excluded keys from both calls and workflows" in {
    // "wf_hello_input" and "cromwell-workflow-id" check against overfiltering: keys should only be filtered relative to
    // workflow and call roots and not relative to arbitrary parts of the JSON object.
    val actual = excludeJson(helloWorldJson, NonEmptyList.of("outputs", "executionEvents", "wf_hello_input", "cromwell-workflow-id")).get
    val expectedMetadata =
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
        |        "inputs": {
        |          "wf_hello_input": "sub world"
        |        },
        |        "end": "2019-07-22T13:32:39.610-04:00",
        |        "attempt": 1,
        |        "start": "2019-07-22T13:32:24.133-04:00",
        |        "subWorkflowId": "ba4a8cec-889e-4cce-9f89-d7df2a12d8d4"
        |      }
        |    ]
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
        |""".stripMargin

    val expected = parseString(expectedMetadata)
    actual shouldEqual expected
  }

  it should "remove excluded keys in subworkflows" in {
    val actual = excludeJson(treblyNestedSubworkflowJson, NonEmptyList.of("executionEvents", "workflowName")).get
    val expected = parseMetadata("excluded_trebly_nested_subworkflow.json")
    actual shouldEqual expected
  }

  it should "keep includes in workflows, paying no mind to that 'labels' field in 'submittedFiles'" in {
    val actual = includeJson(helloWorldJson, NonEmptyList.of("workflowName", "labels", "outputs")).get
    val expectedMetadata = """
      |{
      |  "workflowName" : "main_workflow",
      |  "calls" : {
      |    "main_workflow.wf_hello" : [
      |      {
      |        "shardIndex" : -1,
      |        "outputs" : {
      |          "salutation" : "Hello sub world!"
      |        },
      |        "attempt" : 1
      |      }
      |    ]
      |  },
      |  "outputs" : {
      |    "main_workflow.main_output" : "Hello sub world!"
      |  },
      |  "id" : "757d0bcc-b636-4658-99d4-9b4b3767f1d1",
      |  "labels" : {
      |    "cromwell-workflow-id" : "cromwell-757d0bcc-b636-4658-99d4-9b4b3767f1d1"
      |  }
      |}
      |""".stripMargin

    val expected = parseString(expectedMetadata)
    actual shouldEqual expected
  }

  it should "keep includes in calls and workflows" in {
    val actual = includeJson(helloWorldJson, NonEmptyList.of("workflowName", "executionStatus", "outputs")).get
    val expectedMetadata = """
                             |{
                             |  "workflowName" : "main_workflow",
                             |  "calls" : {
                             |    "main_workflow.wf_hello" : [
                             |      {
                             |        "executionStatus": "Done",
                             |        "shardIndex" : -1,
                             |        "outputs" : {
                             |          "salutation" : "Hello sub world!"
                             |        },
                             |        "attempt" : 1
                             |      }
                             |    ]
                             |  },
                             |  "outputs" : {
                             |    "main_workflow.main_output" : "Hello sub world!"
                             |  },
                             |  "id" : "757d0bcc-b636-4658-99d4-9b4b3767f1d1"
                             |}
                             |""".stripMargin

    val expected = parseString(expectedMetadata)
    actual shouldEqual expected
  }

  // CARBONITE FIXING it should be easy to write a currently-broken version of this like in the Hello World example above.
  it should "keep includes in calls and workflows and subworkflows" in {
    val actual = includeJson(treblyNestedSubworkflowJson, NonEmptyList.of("workflowName", "executionStatus", "outputs")).get
    val expected = parseMetadata("included_trebly_nested_subworkflows.json")
    actual shouldEqual expected
  }

  it should "keep outputs only" in {
    val actual = outputs(helloWorldJson).get
    val expectedMetadata =
      """
        |{
        |  "outputs": {
        |    "main_workflow.main_output": "Hello sub world!"
        |  },
        |  "id": "757d0bcc-b636-4658-99d4-9b4b3767f1d1"
        |}
        |""".stripMargin
    val expected = parseString(expectedMetadata)
    actual shouldEqual expected
  }

  it should "always return a workflow ID on the query workflow" in {
    val expectedMetadata =
      """
        |{
        |  "id" : "757d0bcc-b636-4658-99d4-9b4b3767f1d1"
        |}
        |""".stripMargin

    val expected = parseString(expectedMetadata)
    val actual = includeJson(helloWorldJson, NonEmptyList.of("foo")).toEither.right.get
    actual shouldEqual expected
  }

  it should "do the right thing with logs" in {
    // this actually has no legit logs since its one call is a subworkflow, but we should at least get the workflow ID back.
    val actual = logs(helloWorldJson).get
    val expectedMetadata =
      """
        |{
        |  "id" : "757d0bcc-b636-4658-99d4-9b4b3767f1d1"
        |}
        |""".stripMargin
    val expected = parseString(expectedMetadata)
    actual shouldEqual expected
  }

  // CARBONITE FIXING a subworkflow version of this like the CarboniteMetadataThawingActorSpec would be nice.
  it should "add labels" in {
    val labels = Map(helloWorldJson.workflowId.get -> Map(("new", "label")))
    val newJson = updateLabels(helloWorldJson, labels).get
    val newLabels = newJson.hcursor.downField("labels").keys.get
    assert(newLabels.size === 2)
  }

  it should "replace subworkflow metadata with subworkflow id" in {
    val actualJson = unexpandSubworkflows(treblyNestedSubworkflowJson).get
    val keys = actualJson.hcursor.downField("calls").downField("wf.wf").downArray.keys
    assert(keys.exists(_.exists(_ == "subWorkflowMetadata")) === false)

    val oldCallsArrayCursor = treblyNestedSubworkflowJson.hcursor.downField("calls").downField("wf.wf").downArray
    val oldSubWorkflowId1 = oldCallsArrayCursor.downField("subWorkflowMetadata").get[String]("id").asInstanceOf[Right[DecodingFailure, String]].value
    val oldSubWorkflowId2 = oldCallsArrayCursor.right.downField("subWorkflowMetadata").get[String]("id").asInstanceOf[Right[DecodingFailure, String]].value

    val newCallsArrayCursor = actualJson.hcursor.downField("calls").downField("wf.wf").downArray
    val newSubWorkflowId1 = newCallsArrayCursor.get[String]("subWorkflowId").asInstanceOf[Right[DecodingFailure, String]].value
    val newSubWorkflowId2 = newCallsArrayCursor.right.get[String]("subWorkflowId").asInstanceOf[Right[DecodingFailure, String]].value

    assert(newSubWorkflowId1 === oldSubWorkflowId1)
    assert(newSubWorkflowId2 === oldSubWorkflowId2)

    val thirdArrayElement = newCallsArrayCursor.right.right
    assert(thirdArrayElement.isInstanceOf[FailedCursor])
  }

  it should "properly extract metadata JSON of each subworkflow by id" in {
    val workflowIds = Seq(
      "ba56c1ab-02e0-45f2-97cf-5f91a9138a31",
      "22c6faba-a95a-4e8d-86c3-9a246d7db19b",
      "9e1a2146-f48a-4c04-a589-d66a50dde39b",
      "7d6fba3d-d8e5-43aa-b580-2ace8675ffbd",
      "210b9e04-5606-4231-9cb5-43355d60197d",
      "0e299e7a-bddc-4367-92e2-3e1a61283ca7",
      "be186a2f-b52c-4c6d-96dd-b9a7f16ac526",
      "6382fcbb-fa69-4a6d-bc0f-871013226ad3",
      "55383be2-9a7a-4004-8623-f1cf5a539433",
      "bc649e17-418d-40f6-a145-5a6a8d0c2c5d",
      "0571a73e-1485-4b28-9320-e87036685d61",
      "d9ce3320-727f-42f5-a946-e11177ebd7dd",
      "ffa835a7-68de-4c9d-a777-89f3f7b286dc",
      "540d2d9b-eccc-4e4f-8478-574e4e48f98d"
    )
    workflowIds.foreach { subworkflowId =>
      val extractedSubworkflowJson = extractSubWorkflowMetadata(subworkflowId, treblyNestedSubworkflowJson)
      extractedSubworkflowJson match {
        case Valid(Some(subWorkflowJson)) =>
          assert(true === subWorkflowJson.isInstanceOf[Json])
          assert(subworkflowId === subWorkflowJson.workflowId.get.toString)
        case Valid(None) => fail(s"Subworkflow not found for id $subworkflowId")
        case Invalid(e) => fail(e.toList.mkString("\n"))
      }
    }
  }

  it should "start including and excluding keys from both (but only!) workflow- and call-roots" in {
    val newJson = includeJson(helloWorldJson, NonEmptyList.of("start")).get
    val expectedJsonRaw =
      """
        |{
        |  "calls": {
        |    "main_workflow.wf_hello": [
        |      {
        |        "shardIndex": -1,
        |        "attempt": 1,
        |        "start": "2019-07-22T13:32:24.133-04:00"
        |      }
        |    ]
        |  },
        |  "id": "757d0bcc-b636-4658-99d4-9b4b3767f1d1",
        |  "start": "2019-07-22T13:32:20.434-04:00"
        |}
      """.stripMargin
    val expectedJson = parse(expectedJsonRaw)

    newJson shouldEqual expectedJson.right.get
  }

  it should "support a Job Manageresque query on a simple workflow" in {
    import JobManagerKeys._

    val actual = includeExcludeJson(helloWorldJson, includeKeys, excludeKeys).toEither.right.get
    val expected = parseMetadata("jmui_hello_world.json")
    actual shouldEqual expected
  }

  it should "support a Job Manageresque query on a workflow having trebly nested subworkflows" in {
    import JobManagerKeys._

    val actual = includeExcludeJson(parseMetadata("trebly_nested_subworkflow_papiv2_cached.json"), includeKeys, excludeKeys).toEither.right.get
    val expected = parseMetadata("trebly_nested_jmuiesque.json")
    actual shouldEqual expected
  }

  it should "handle includes nested multiple levels" in {
    val includeKeys = NonEmptyList.of("callCaching:hashes:backend name", "callCaching:hashes:input:String addressee")
    val actual = includeJson(helloWorldPapiV2, includeKeys).get
    val expectedJson =
      """
        |{
        |  "calls" : {
        |    "wf_hello.hello" : [
        |      {
        |        "callCaching" : {
        |          "hashes": {
        |            "input" : {
        |              "String addressee" : "FF3B11B1B0C96178D517D6EAF44F5519"
        |            },
        |            "backend name" : "36EF4A8AB268D1A1C74D8108C93D48ED"
        |          }
        |        },
        |        "shardIndex" : -1,
        |        "attempt" : 1
        |      }
        |    ]
        |  },
        |  "id" : "d53a063a-e8b7-403f-a400-a85f089a8928"
        |}
        |""".stripMargin

    val expected = parseString(expectedJson)
    actual shouldEqual expected
  }

  it should "handle includes specifying nonexistent keys" in {
    val includeKeys = NonEmptyList.of("callCaching:hashes:nonExistent", "callCaching:phonyKey")
    val actual = includeJson(helloWorldPapiV2, includeKeys).get
    val expectedJson =
      """
        |{
        |  "id": "d53a063a-e8b7-403f-a400-a85f089a8928"
        |}
        |""".stripMargin
    val expected = parseString(expectedJson)
    actual shouldEqual expected
  }

  it should "handle excludes specifying nonexistent keys" in {
    val includeKeys = NonEmptyList.of("backendLabels")
    val excludeKeys = NonEmptyList.of("backendLabels:foo")

    val actual = includeExcludeJson(helloWorldPapiV2, includeKeys, excludeKeys).get
    val expectedJson =
      """
        |{
        |  "id": "d53a063a-e8b7-403f-a400-a85f089a8928",
        |  "calls" : {
        |    "wf_hello.hello" : [
        |      {
        |        "backendLabels": {
        |          "cromwell-workflow-id": "cromwell-d53a063a-e8b7-403f-a400-a85f089a8928",
        |          "wdl-task-name": "hello"
        |        },
        |        "shardIndex" : -1,
        |        "attempt" : 1
        |      }
        |    ]
        |  }
        |}
        |""".stripMargin
    val expected = parseString(expectedJson)
    actual shouldEqual expected
  }

  it should "be able to deal with some nested include / exclude stuff" in {
    val includeKeys = NonEmptyList.of("callCaching")
    val excludeKeys = NonEmptyList.of("callCaching:hashes:input:String addressee", "callCaching:hashes:runtime attribute:docker")

    val actual = includeExcludeJson(helloWorldPapiV2, includeKeys, excludeKeys).get
    val expectedJson =
      """
        |{
        |  "calls": {
        |    "wf_hello.hello": [
        |      {
        |        "callCaching": {
        |          "allowResultReuse": true,
        |          "effectiveCallCachingMode": "ReadAndWriteCache",
        |          "hashes": {
        |            "backend name": "36EF4A8AB268D1A1C74D8108C93D48ED",
        |            "command template": "4EAADE3CD5D558C5A6CFA4FD101A1486",
        |            "input count": "C4CA4238A0B923820DCC509A6F75849B",
        |            "output count": "C4CA4238A0B923820DCC509A6F75849B",
        |            "output expression": {
        |              "String salutation": "0183144CF6617D5341681C6B2F756046"
        |            },
        |            "runtime attribute": {
        |              "continueOnReturnCode": "CFCD208495D565EF66E7DFF9F98764DA",
        |              "failOnStderr": "68934A3E9455FA72420237EB05902327"
        |            }
        |          },
        |          "hit": false,
        |          "result": "Cache Miss"
        |        },
        |        "attempt": 1,
        |        "shardIndex": -1
        |      }
        |    ]
        |  },
        |  "id": "d53a063a-e8b7-403f-a400-a85f089a8928"
        |}
        |""".stripMargin

    val expected = parseString(expectedJson)
    actual shouldEqual expected
  }

  it should "exclude keys which are listed to both 'include' and 'exclude'" in {
    val includeKeys = NonEmptyList.of("callCaching", "backendStatus")
    val excludeKeys = NonEmptyList.of("callCaching", "callCaching:hashes:runtime attribute:docker")

    val actual = includeExcludeJson(helloWorldPapiV2, includeKeys, excludeKeys).get
    val expectedJson =
      """
        |{
        |  "calls": {
        |    "wf_hello.hello": [
        |      {
        |        "backendStatus": "Success",
        |        "attempt": 1,
        |        "shardIndex": -1
        |      }
        |    ]
        |  },
        |  "id": "d53a063a-e8b7-403f-a400-a85f089a8928"
        |}
        |""".stripMargin

    val expected = parseString(expectedJson)
    actual shouldEqual expected
  }

  it should "remove excluded sections even if a sub-section was included" in {
    val includeKeys = NonEmptyList.of("callCaching:hashes", "backendStatus")
    val excludeKeys = NonEmptyList.of("callCaching")

    val actual = includeExcludeJson(helloWorldPapiV2, includeKeys, excludeKeys).get
    val expectedJson =
      """
        |{
        |  "calls": {
        |    "wf_hello.hello": [
        |      {
        |        "backendStatus": "Success",
        |        "attempt": 1,
        |        "shardIndex": -1
        |      }
        |    ]
        |  },
        |  "id": "d53a063a-e8b7-403f-a400-a85f089a8928"
        |}
        |""".stripMargin

    val expected = parseString(expectedJson)
    actual shouldEqual expected
  }

  it should "ignore attempts to exclude ids, attempts, and shard indices" in {
    val includeKeys = None
    val excludeKeys = Option(NonEmptyList.of("id", "attempt", "shardIndex"))

    val actual = includeExcludeJson(helloWorldPapiV2, includeKeys, excludeKeys).get
    val expected = helloWorldPapiV2

    actual shouldEqual expected
  }

  it should "be able to exclude calls, if so asked" in {
    val includeKeys = None
    val excludeKeys = Option(NonEmptyList.of("calls"))

    val actual = includeExcludeJson(helloWorldPapiV2, includeKeys, excludeKeys).get
    val expected = helloWorldPapiV2.mapObject(_.remove("calls"))
    actual shouldEqual expected
  }
}

object JsonEditorSpec {

  private def parseString(string: String): Json = parse(string).right.get

  private def parseMetadata(filename: String): Json = {
    parseString(
      Source
        .fromInputStream(Thread
          .currentThread
          .getContextClassLoader
          .getResourceAsStream(filename)
        ).mkString)
  }

  val treblyNestedSubworkflowJson: Json = parseMetadata("trebly_nested_subworkflows.json")
  val treblyNestedSubworkflowCachedJson: Json = parseMetadata("trebly_nested_subworkflow_papiv2_cached.json")
  val helloWorldJson: Json = parseMetadata("hello_world.json")
  val helloWorldPapiV2: Json = parseMetadata("hello_papiv2.json")

  implicit class EnhancedErrorOr[A](val e: ErrorOr[A]) extends AnyVal {
    def get: A = e.toEither.right.get
  }

  object JobManagerKeys {
    val includeKeys = NonEmptyList.of(
      "attempt",
      "backendLogs:log",
      "callCaching:hit",
      "callRoot",
      "end",
      "executionStatus",
      "failures",
      "inputs",
      "jobId",
      "outputs",
      "shardIndex",
      "start",
      "calls",
      "description",
      "executionEvents",
      "labels",
      "parentWorkflowId",
      "returnCode",
      "status",
      "submission",
      "subWorkflowId",
      "workflowName"
    )

    val excludeKeys = NonEmptyList.of("callCaching:hitFailures")
  }
}
