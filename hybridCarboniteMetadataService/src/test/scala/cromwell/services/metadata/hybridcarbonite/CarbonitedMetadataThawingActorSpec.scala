package cromwell.services.metadata.hybridcarbonite

import java.util.UUID

import akka.stream.ActorMaterializer
import akka.testkit.{TestActorRef, TestProbe}
import com.typesafe.config.ConfigFactory
import common.validation.Validation._
import cromwell.core.io.IoContentAsStringCommand
import cromwell.core.io.IoPromiseProxyActor.IoCommandWithPromise
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.services.metadata.MetadataService.{GetLabels, LabelLookupResponse}
import cromwell.services.metadata.hybridcarbonite.CarbonitedMetadataThawingActor.{ThawCarboniteFailed, ThawCarboniteSucceeded, ThawCarbonitedMetadata}
import cromwell.services.metadata.hybridcarbonite.CarbonitedMetadataThawingActorSpec._
import org.scalatest.{FlatSpecLike, Matchers}
import io.circe.parser._

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

class CarbonitedMetadataThawingActorSpec extends TestKitSuite("CarbonitedMetadataThawingActorSpec") with FlatSpecLike with Matchers {

  implicit val ec: ExecutionContext = system.dispatcher

  val carboniterConfig = HybridCarboniteConfig.make(ConfigFactory.parseString(
    """bucket = "carbonite-test-bucket"
      |filesystems {
      |  gcs {
      |    # A reference to the auth to use for storing and retrieving metadata:
      |    auth = "application-default"
      |  }
      |}""".stripMargin)).unsafe("Make config file")

  val serviceRegistryActor = TestProbe()
  val ioActor = TestProbe()

  implicit val materializer: ActorMaterializer = ActorMaterializer()

  it should "receive a message from GCS" in {

    val clientProbe = TestProbe()

    val actorUnderTest = TestActorRef(new CarbonitedMetadataThawingActor(carboniterConfig, serviceRegistryActor.ref, ioActor.ref), "ThawingActor")

    clientProbe.send(actorUnderTest, ThawCarbonitedMetadata(workflowId))

    serviceRegistryActor.expectMsg(GetLabels(workflowId))
    serviceRegistryActor.send(actorUnderTest, LabelLookupResponse(workflowId, Map("bob loblaw" -> "law blog")))

    ioActor.expectMsgPF(max = 5.seconds) {
      case command @ IoCommandWithPromise(iocasc: IoContentAsStringCommand, _) if iocasc.file.pathAsString.contains(workflowId.toString) =>
        command.promise.success(rawMetadataSample)
    }

    clientProbe.expectMsgPF(max = 5.seconds) {
      case ThawCarboniteSucceeded(str) => parse(str) should be(parse(augmentedMetadataSample))
      case ThawCarboniteFailed(reason) => fail(reason)
    }
  }
}

object CarbonitedMetadataThawingActorSpec {
  val workflowId = WorkflowId(UUID.fromString("2ce544a0-4c0d-4cc9-8a0b-b412bb1e5f82"))

  val rawMetadataSample = sampleMetadataContent("")
  val augmentedMetadataSample = sampleMetadataContent("""
                                          |"labels" : {
                                          |    "bob loblaw" : "law blog"
                                          |},
                                          |  """.stripMargin)

  def sampleMetadataContent(labelsContent: String) = s"""{$labelsContent
                                                        |  "workflowName": "helloWorldWf",
                                                        |  "workflowProcessingEvents": [
                                                        |    {
                                                        |      "cromwellId": "cromid-78ea393",
                                                        |      "description": "PickedUp",
                                                        |      "timestamp": "2019-07-31T20:19:02.165Z",
                                                        |      "cromwellVersion": "43-e9c0f4c-SNAP"
                                                        |    },
                                                        |    {
                                                        |      "cromwellId": "cromid-78ea393",
                                                        |      "description": "Finished",
                                                        |      "timestamp": "2019-07-31T20:19:35.058Z",
                                                        |      "cromwellVersion": "43-e9c0f4c-SNAP"
                                                        |    }
                                                        |  ],
                                                        |  "actualWorkflowLanguageVersion": "draft-2",
                                                        |  "submittedFiles": {
                                                        |    "workflow": "task helloWorld {\\n    \\n    command { echo \\"hello, world\\" }\\n\\n    output { String s = read_string(stdout()) }\\n    \\n    runtime {\\n      docker: \\"ubuntu:latest\\"\\n    }\\n}\\n\\nworkflow helloWorldWf {\\n    call helloWorld\\n}\\n",
                                                        |    "root": "",
                                                        |    "options": "{\\n\\n}",
                                                        |    "inputs": "{}",
                                                        |    "workflowUrl": "",
                                                        |    "labels": "{}"
                                                        |  },
                                                        |  "calls": {
                                                        |    "helloWorldWf.helloWorld": [
                                                        |      {
                                                        |        "executionStatus": "Done",
                                                        |        "stdout": "gs://execution-bucket/$workflowId/call-name/stdout",
                                                        |        "backendStatus": "Done",
                                                        |        "compressedDockerSize": 26723061,
                                                        |        "commandLine": "echo \\"hello, world\\"",
                                                        |        "shardIndex": -1,
                                                        |        "outputs": {
                                                        |          "s": "hello, world"
                                                        |        },
                                                        |        "runtimeAttributes": {
                                                        |          "docker": "ubuntu:latest",
                                                        |          "failOnStderr": "false",
                                                        |          "maxRetries": "0",
                                                        |          "continueOnReturnCode": "0"
                                                        |        },
                                                        |        "callCaching": {
                                                        |          "allowResultReuse": true,
                                                        |          "hit": false,
                                                        |          "result": "Cache Miss",
                                                        |          "hashes": {
                                                        |            "output count": "C4CA4238A0B923820DCC509A6F75849B",
                                                        |            "runtime attribute": {
                                                        |              "docker": "CFE62D82138579B4B9F4EE81EFC5745E",
                                                        |              "continueOnReturnCode": "CFCD208495D565EF66E7DFF9F98764DA",
                                                        |              "failOnStderr": "68934A3E9455FA72420237EB05902327"
                                                        |            },
                                                        |            "output expression": {
                                                        |              "String s": "0183144CF6617D5341681C6B2F756046"
                                                        |            },
                                                        |            "input count": "CFCD208495D565EF66E7DFF9F98764DA",
                                                        |            "backend name": "509820290D57F333403F490DDE7316F4",
                                                        |            "command template": "720B86AE4F6DFD9C11FD9E01AEC6FDF9"
                                                        |          },
                                                        |          "effectiveCallCachingMode": "ReadAndWriteCache"
                                                        |        },
                                                        |        "inputs": {},
                                                        |        "returnCode": 0,
                                                        |        "jobId": "5399",
                                                        |        "backend": "Local",
                                                        |        "end": "2019-07-31T20:19:33.251Z",
                                                        |        "dockerImageUsed": "ubuntu@sha256:c303f19cfe9ee92badbbbd7567bc1ca47789f79303ddcef56f77687d4744cd7a",
                                                        |        "stderr": "gs://execution-bucket/$workflowId/call-name/stderr",
                                                        |        "callRoot": "gs://execution-bucket/$workflowId/call-name",
                                                        |        "attempt": 1,
                                                        |        "executionEvents": [
                                                        |          {
                                                        |            "startTime": "2019-07-31T20:19:05.344Z",
                                                        |            "description": "PreparingJob",
                                                        |            "endTime": "2019-07-31T20:19:06.450Z"
                                                        |          },
                                                        |          {
                                                        |            "startTime": "2019-07-31T20:19:06.532Z",
                                                        |            "description": "RunningJob",
                                                        |            "endTime": "2019-07-31T20:19:30.217Z"
                                                        |          },
                                                        |          {
                                                        |            "startTime": "2019-07-31T20:19:04.447Z",
                                                        |            "description": "RequestingExecutionToken",
                                                        |            "endTime": "2019-07-31T20:19:05.328Z"
                                                        |          },
                                                        |          {
                                                        |            "startTime": "2019-07-31T20:19:04.416Z",
                                                        |            "description": "Pending",
                                                        |            "endTime": "2019-07-31T20:19:04.447Z"
                                                        |          },
                                                        |          {
                                                        |            "startTime": "2019-07-31T20:19:30.217Z",
                                                        |            "description": "UpdatingCallCache",
                                                        |            "endTime": "2019-07-31T20:19:32.290Z"
                                                        |          },
                                                        |          {
                                                        |            "startTime": "2019-07-31T20:19:05.328Z",
                                                        |            "description": "WaitingForValueStore",
                                                        |            "endTime": "2019-07-31T20:19:05.344Z"
                                                        |          },
                                                        |          {
                                                        |            "startTime": "2019-07-31T20:19:32.290Z",
                                                        |            "description": "UpdatingJobStore",
                                                        |            "endTime": "2019-07-31T20:19:33.251Z"
                                                        |          },
                                                        |          {
                                                        |            "startTime": "2019-07-31T20:19:06.450Z",
                                                        |            "description": "CheckingCallCache",
                                                        |            "endTime": "2019-07-31T20:19:06.532Z"
                                                        |          }
                                                        |        ],
                                                        |        "start": "2019-07-31T20:19:04.392Z"
                                                        |      }
                                                        |    ]
                                                        |  },
                                                        |  "outputs": {
                                                        |    "helloWorldWf.helloWorld.s": "hello, world"
                                                        |  },
                                                        |  "workflowRoot": "gs://execution-bucket/$workflowId",
                                                        |  "actualWorkflowLanguage": "WDL",
                                                        |  "id": "$workflowId",
                                                        |  "inputs": {},
                                                        |  "submission": "2019-07-31T20:19:02.067Z",
                                                        |  "status": "Succeeded",
                                                        |  "end": "2019-07-31T20:19:35.057Z",
                                                        |  "start": "2019-07-31T20:19:02.228Z"
                                                        |}""".stripMargin
}
