package centaur.json

import centaur.json.JsonUtils.EnhancedJsValue
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json._

class JsonUtilsSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "JsonUtils"

  it should "correctly flatten simple sample metadata" in {
    val metadata =
      """{
        | "actualWorkflowLanguage": "WDL",
        | "actualWorkflowLanguageVersion": "draft-2",
        | "id": "5abfaa90-570f-48d4-a35b-81d5ad4ea0f7",
        | "status": "Succeeded",
        | "workflowName": "wf_hello",
        | "inputs": {
        |    "wf_hello.hello.addressee": "Mr. Bean"
        |  },
        |  "outputs": {
        |    "wf_hello.hello.salutation": "Hello Mr. Bean!"
        |  }
        |}""".stripMargin

    val expectedFlattenedMetadata: Map[String, JsValue] = Map(
      "actualWorkflowLanguage" -> "WDL",
      "actualWorkflowLanguageVersion" -> "draft-2",
      "id" -> "5abfaa90-570f-48d4-a35b-81d5ad4ea0f7",
      "status" -> "Succeeded",
      "workflowName" -> "wf_hello",
      "inputs.wf_hello.hello.addressee" -> "Mr. Bean",
      "outputs.wf_hello.hello.salutation" -> "Hello Mr. Bean!"
    ).map(x => (x._1, JsString(x._2)))

    val actualFlattenedMetadata: Map[String, JsValue] = metadata.parseJson.asJsObject.flatten().fields

    actualFlattenedMetadata should contain theSameElementsAs expectedFlattenedMetadata
  }

  it should "correctly flatten call metadata with single task with 1 attempt" in {
    val metadata =
      """{
        | "id": "5abfaa90-570f-48d4-a35b-81d5ad4ea0f7",
        | "status": "Succeeded",
        | "workflowName": "wf_hello",
        | "calls": {
        |    "wf_hello.hello": [
        |      {
        |        "attempt": 1,
        |        "backend": "Papiv2",
        |        "backendStatus": "Success",
        |        "executionEvents": [
        |          {
        |            "description": "step 1",
        |            "endTime": "2021-01-21T21:46:39.748Z",
        |            "startTime": "2021-01-21T21:46:20.456Z"
        |          },
        |          {
        |            "description": "step 2",
        |            "endTime": "2021-01-21T21:46:41.685Z",
        |            "startTime": "2021-01-21T21:46:39.748Z"
        |          },
        |          {
        |            "description": "step 3",
        |            "endTime": "2021-01-21T21:47:51.946Z",
        |            "startTime": "2021-01-21T21:47:51.771Z"
        |          }
        |        ],
        |        "executionStatus": "Done",
        |        "inputs": {
        |          "addressee": "Mr. Bean"
        |        },
        |        "outputs": {
        |          "salutation": "Hello Mr. Bean!"
        |        },
        |        "runtimeAttributes": {
        |          "bootDiskSizeGb": "10",
        |          "continueOnReturnCode": "0",
        |          "maxRetries": "0"
        |        },
        |        "shardIndex": -1
        |      }
        |    ]
        |  }
        |}""".stripMargin

    val expectedFlattenedMetadata: Map[String, JsValue] = Map(
      "calls.wf_hello.hello.attempt" -> JsNumber(1),
      "calls.wf_hello.hello.shardIndex" -> JsNumber(-1)
    ) ++ Map(
      "id" -> "5abfaa90-570f-48d4-a35b-81d5ad4ea0f7",
      "status" -> "Succeeded",
      "workflowName" -> "wf_hello",
      "calls.wf_hello.hello.backend" -> "Papiv2",
      "calls.wf_hello.hello.backendStatus" -> "Success",
      "calls.wf_hello.hello.executionEvents.0.description" -> "step 1",
      "calls.wf_hello.hello.executionEvents.0.endTime" -> "2021-01-21T21:46:39.748Z",
      "calls.wf_hello.hello.executionEvents.0.startTime" -> "2021-01-21T21:46:20.456Z",
      "calls.wf_hello.hello.executionEvents.1.description" -> "step 2",
      "calls.wf_hello.hello.executionEvents.1.endTime" -> "2021-01-21T21:46:41.685Z",
      "calls.wf_hello.hello.executionEvents.1.startTime" -> "2021-01-21T21:46:39.748Z",
      "calls.wf_hello.hello.executionEvents.2.description" -> "step 3",
      "calls.wf_hello.hello.executionEvents.2.endTime" -> "2021-01-21T21:47:51.946Z",
      "calls.wf_hello.hello.executionEvents.2.startTime" -> "2021-01-21T21:47:51.771Z",
      "calls.wf_hello.hello.executionStatus" -> "Done",
      "calls.wf_hello.hello.inputs.addressee" -> "Mr. Bean",
      "calls.wf_hello.hello.outputs.salutation" -> "Hello Mr. Bean!",
      "calls.wf_hello.hello.runtimeAttributes.bootDiskSizeGb" -> "10",
      "calls.wf_hello.hello.runtimeAttributes.continueOnReturnCode" -> "0",
      "calls.wf_hello.hello.runtimeAttributes.maxRetries" -> "0",
    ).map(x => (x._1, JsString(x._2)))

    val actualFlattenedMetadata: Map[String, JsValue] = metadata.parseJson.asJsObject.flatten().fields

    actualFlattenedMetadata should contain theSameElementsAs expectedFlattenedMetadata
  }

  it should "correctly flatten call metadata with 2 tasks with 1 attempt each" in {
    val metadata =
      """{
        | "workflowName": "wf_hello",
        | "calls": {
        |    "wf_hello.task1": [
        |      {
        |        "attempt": 1,
        |        "executionEvents": [
        |          {
        |            "description": "task 1 step 1"
        |          }
        |        ],
        |        "executionStatus": "Done",
        |        "runtimeAttributes": {
        |          "bootDiskSizeGb": "10"
        |        },
        |        "shardIndex": -1
        |      }
        |    ],
        |    "wf_hello.task2": [
        |      {
        |        "attempt": 1,
        |        "executionEvents": [
        |          {
        |            "description": "task 2 step 1"
        |          }
        |        ],
        |        "executionStatus": "Failed",
        |        "runtimeAttributes": {
        |          "bootDiskSizeGb": "10"
        |        },
        |        "shardIndex": -1
        |      }
        |    ]
        |  }
        |}""".stripMargin

    val expectedFlattenedMetadata: Map[String, JsValue] = Map(
      "calls.wf_hello.task1.attempt" -> JsNumber(1),
      "calls.wf_hello.task2.attempt" -> JsNumber(1),
      "calls.wf_hello.task1.shardIndex" -> JsNumber(-1),
      "calls.wf_hello.task2.shardIndex" -> JsNumber(-1)
    ) ++ Map(
      "workflowName" -> "wf_hello",
      "calls.wf_hello.task1.executionStatus" -> "Done",
      "calls.wf_hello.task2.executionStatus" -> "Failed",
      "calls.wf_hello.task1.executionEvents.0.description" -> "task 1 step 1",
      "calls.wf_hello.task2.executionEvents.0.description" -> "task 2 step 1",
      "calls.wf_hello.task1.runtimeAttributes.bootDiskSizeGb" -> "10",
      "calls.wf_hello.task2.runtimeAttributes.bootDiskSizeGb" -> "10",
    ).map(x => (x._1, JsString(x._2)))

    val actualFlattenedMetadata: Map[String, JsValue] = metadata.parseJson.asJsObject.flatten().fields

    actualFlattenedMetadata should contain theSameElementsAs expectedFlattenedMetadata
  }

  it should "correctly flatten call metadata with single task, multiple shards and multiple attempts each" in {
    val metadata =
      """{
        | "id": "5abfaa90-570f-48d4-a35b-81d5ad4ea0f7",
        | "status": "Succeeded",
        | "workflowName": "wf_hello",
        | "calls": {
        |    "wf_hello.hello": [
        |      {
        |        "attempt": 1,
        |        "executionEvents": [
        |          {
        |            "description": "shard 0 attempt 1 step 1"
        |          }
        |        ],
        |        "runtimeAttributes": {
        |          "memory": "1 GB"
        |        },
        |        "shardIndex": 0
        |      },
        |      {
        |        "attempt": 2,
        |        "executionEvents": [
        |          {
        |            "description": "shard 0 attempt 2 step 1"
        |          }
        |        ],
        |        "runtimeAttributes": {
        |          "memory": "1.1 GB"
        |        },
        |        "shardIndex": 0
        |      },
        |      {
        |        "attempt": 1,
        |        "executionEvents": [
        |          {
        |            "description": "shard 1 attempt 1 step 1"
        |          }
        |        ],
        |        "runtimeAttributes": {
        |          "memory": "1 GB"
        |        },
        |        "shardIndex": 1
        |      },
        |      {
        |        "attempt": 2,
        |        "executionEvents": [
        |          {
        |            "description": "shard 1 attempt 2 step 1"
        |          }
        |        ],
        |        "runtimeAttributes": {
        |          "memory": "1.1 GB"
        |        },
        |        "shardIndex": 1
        |      }
        |    ]
        |  }
        |}""".stripMargin

    val expectedFlattenedMetadata: Map[String, JsValue] = Map(
      "wf_hello.hello.0.attempt" -> 2,
      "wf_hello.hello.0.1.attempt" -> 1,
      "wf_hello.hello.0.2.attempt" -> 2,
      "wf_hello.hello.1.attempt" -> 2,
      "wf_hello.hello.1.1.attempt" -> 1,
      "wf_hello.hello.1.2.attempt" -> 2,
      "wf_hello.hello.0.shardIndex" -> 0,
      "wf_hello.hello.0.1.shardIndex" -> 0,
      "wf_hello.hello.0.2.shardIndex" -> 0,
      "wf_hello.hello.1.shardIndex" -> 1,
      "wf_hello.hello.1.1.shardIndex" -> 1,
      "wf_hello.hello.1.2.shardIndex" -> 1,
    ).map(x => (x._1, JsNumber(x._2))) ++ Map(
      "id" -> "5abfaa90-570f-48d4-a35b-81d5ad4ea0f7",
      "status" -> "Succeeded",
      "workflowName" -> "wf_hello",
      "wf_hello.hello.0.executionEvents.0.description" -> "shard 0 attempt 2 step 1",
      "wf_hello.hello.0.1.executionEvents.0.description" -> "shard 0 attempt 1 step 1",
      "wf_hello.hello.0.2.executionEvents.0.description" -> "shard 0 attempt 2 step 1",
      "wf_hello.hello.1.executionEvents.0.description" -> "shard 1 attempt 2 step 1",
      "wf_hello.hello.1.1.executionEvents.0.description" -> "shard 1 attempt 1 step 1",
      "wf_hello.hello.1.2.executionEvents.0.description" -> "shard 1 attempt 2 step 1",
      "wf_hello.hello.0.runtimeAttributes.memory" -> "1.1 GB",
      "wf_hello.hello.0.1.runtimeAttributes.memory" -> "1 GB",
      "wf_hello.hello.0.2.runtimeAttributes.memory" -> "1.1 GB",
      "wf_hello.hello.1.runtimeAttributes.memory" -> "1.1 GB",
      "wf_hello.hello.1.1.runtimeAttributes.memory" -> "1 GB",
      "wf_hello.hello.1.2.runtimeAttributes.memory" -> "1.1 GB",
    ).map(x => (x._1, JsString(x._2)))

    val actualFlattenedMetadata: Map[String, JsValue] = metadata.parseJson.asJsObject.flatten().fields

    actualFlattenedMetadata should contain theSameElementsAs expectedFlattenedMetadata
  }
}
