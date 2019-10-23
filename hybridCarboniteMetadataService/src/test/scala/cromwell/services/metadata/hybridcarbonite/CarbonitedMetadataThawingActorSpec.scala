package cromwell.services.metadata.hybridcarbonite

import akka.stream.ActorMaterializer
import akka.testkit.{TestActorRef, TestProbe}
import com.fasterxml.jackson.databind.{JsonNode, ObjectMapper}
import com.typesafe.config.ConfigFactory
import common.validation.Validation._
import cromwell.core.io.IoContentAsStringCommand
import cromwell.core.io.IoPromiseProxyActor.IoCommandWithPromise
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.services.metadata.MetadataQuery
import cromwell.services.metadata.MetadataService.{GetMetadataAction, GetRootAndSubworkflowLabels, RootAndSubworkflowLabelsLookupResponse}
import cromwell.services.metadata.hybridcarbonite.CarbonitedMetadataThawingActorSpec.{workflowId, _}
import cromwell.services.{BuiltMetadataResponse, FailedMetadataResponse}
import io.circe.parser._
import net.thisptr.jackson.jq._
import org.scalatest.{FlatSpecLike, Matchers}

import scala.collection.JavaConverters._
import scala.collection.mutable.ArrayBuffer
import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.io.Source
import scala.language.postfixOps


class CarbonitedMetadataThawingActorSpec extends TestKitSuite("CarbonitedMetadataThawingActorSpec") with FlatSpecLike with Matchers {

  implicit val ec: ExecutionContext = system.dispatcher

  val carboniterConfig = HybridCarboniteConfig.parseConfig(ConfigFactory.parseString(
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
    clientProbe.send(actorUnderTest, GetMetadataAction(
      MetadataQuery(
        workflowId = workflowId,
        jobKey = None,
        key = None,
        includeKeysOption = None,
        excludeKeysOption = None,
        expandSubWorkflows = true
      )))

    serviceRegistryActor.expectMsg(GetRootAndSubworkflowLabels(workflowId))
    serviceRegistryActor.send(actorUnderTest, RootAndSubworkflowLabelsLookupResponse(workflowId, Map(WorkflowId(workflowId.id) -> Map("bob loblaw" -> "law blog"))))

    ioActor.expectMsgPF(max = 5.seconds) {
      case command @ IoCommandWithPromise(iocasc: IoContentAsStringCommand, _) if iocasc.file.pathAsString.contains(workflowId.toString) =>
        command.promise.success(rawMetadataSample)
    }

    clientProbe.expectMsgPF(max = 5.seconds) {
      case BuiltMetadataResponse(_, jsObject) => parse(jsObject.compactPrint) should be(parse(augmentedMetadataSample))
      case FailedMetadataResponse(_, reason) => fail(reason)
    }
  }

  it should "add labels to subworkflows" in {

    val rootWorkflowId = WorkflowId.fromString("e0c918bb-695e-458a-ab37-be33db7bf721")
    // The subworkflow IDs in the order they appear in the metadata JSON.
    val subWorkflowIds = List(
      "be186a2f-b52c-4c6d-96dd-b9a7f16ac526",
      "9e1a2146-f48a-4c04-a589-d66a50dde39b",
      "ba56c1ab-02e0-45f2-97cf-5f91a9138a31",
      "22c6faba-a95a-4e8d-86c3-9a246d7db19b",
      "0e299e7a-bddc-4367-92e2-3e1a61283ca7",
      "7d6fba3d-d8e5-43aa-b580-2ace8675ffbd",
      "210b9e04-5606-4231-9cb5-43355d60197d",
      "540d2d9b-eccc-4e4f-8478-574e4e48f98d",
      "bc649e17-418d-40f6-a145-5a6a8d0c2c5d",
      "6382fcbb-fa69-4a6d-bc0f-871013226ad3",
      "55383be2-9a7a-4004-8623-f1cf5a539433",
      "ffa835a7-68de-4c9d-a777-89f3f7b286dc",
      "0571a73e-1485-4b28-9320-e87036685d61",
      "d9ce3320-727f-42f5-a946-e11177ebd7dd") map WorkflowId.fromString

    val labelsForUpdate: Map[WorkflowId, Map[String, String]] = (rootWorkflowId :: subWorkflowIds) map { id => id -> Map("short_id" -> id.shortString) } toMap

    val action = GetMetadataAction(
      MetadataQuery(
        workflowId = rootWorkflowId,
        jobKey = None,
        key = None,
        includeKeysOption = None,
        excludeKeysOption = None,
        expandSubWorkflows = true
      ))

    val clientProbe = TestProbe()
    val actorUnderTest = TestActorRef(new CarbonitedMetadataThawingActor(carboniterConfig, serviceRegistryActor.ref, ioActor.ref), "ThawingActor")

    clientProbe.send(actorUnderTest, action)
    serviceRegistryActor.expectMsg(GetRootAndSubworkflowLabels(rootWorkflowId))
    serviceRegistryActor.send(actorUnderTest, RootAndSubworkflowLabelsLookupResponse(rootWorkflowId, labelsForUpdate))

    val metadataWithSubworkflows = Source.fromInputStream(Thread.currentThread.getContextClassLoader.getResourceAsStream("metadata_with_subworkflows.json")).mkString

    val scope = Scope.newEmptyScope()
    BuiltinFunctionLoader.getInstance.loadFunctions(Versions.JQ_1_5, scope)
    val objectMapper = new ObjectMapper()

    def outputBuilder(buf: ArrayBuffer[JsonNode]): Output = out => buf.append(out)

    val rootWorkflowJsonNodesBeforeUpdate = new ArrayBuffer[JsonNode]()
    val beforeUpdateOutput = outputBuilder(rootWorkflowJsonNodesBeforeUpdate)

    val rootWorkflowLabelsQuery = JsonQuery.compile(".labels", Versions.JQ_1_5)
    val rootWorkflowJsonNodeBeforeUpdate = objectMapper.readTree(metadataWithSubworkflows)
    rootWorkflowLabelsQuery.apply(scope, rootWorkflowJsonNodeBeforeUpdate, beforeUpdateOutput)
    rootWorkflowJsonNodesBeforeUpdate.size shouldBe 1
    rootWorkflowJsonNodesBeforeUpdate.head.size() shouldBe 1
    val actualRootLabelsBeforeUpdate = rootWorkflowJsonNodesBeforeUpdate.head.fields().asScala.toList map { e => e.getKey -> e.getValue.textValue() } toMap
    val expectedRootLabelsBeforeUpdate = Map("cromwell-workflow-id" -> ("cromwell-" + rootWorkflowId.toString))
    actualRootLabelsBeforeUpdate shouldEqual expectedRootLabelsBeforeUpdate

    ioActor.expectMsgPF(max = 5.seconds) {
      case command @ IoCommandWithPromise(iocasc: IoContentAsStringCommand, _) if iocasc.file.pathAsString.contains(rootWorkflowId.toString) =>
        command.promise.success(metadataWithSubworkflows)
    }

    clientProbe.expectMsgPF(max = 5.seconds) {
      case BuiltMetadataResponse(_, metadataAfterThawing) =>
        val rootNodesAfterUpdate = new ArrayBuffer[JsonNode]()
        val rootOutputAfterUpdate = outputBuilder(rootNodesAfterUpdate)

        val jsonNodeAfterUpdate = objectMapper.readTree(metadataAfterThawing.compactPrint)
        rootWorkflowLabelsQuery.apply(scope, jsonNodeAfterUpdate, rootOutputAfterUpdate)
        val actualRootLabelsAfterUpdate = rootNodesAfterUpdate.head.fields().asScala.toList map { e => e.getKey -> e.getValue.textValue() } toMap
        val expectedRootLabelsAfterUpdate = Map(
          "cromwell-workflow-id" -> ("cromwell-" + rootWorkflowId.toString),
          "short_id" -> rootWorkflowId.shortString
        )

        actualRootLabelsAfterUpdate shouldEqual expectedRootLabelsAfterUpdate

        val subworkflowNodesAfterUpdate = new ArrayBuffer[JsonNode]()
        val subworkflowOutputAfterUpdate = outputBuilder(subworkflowNodesAfterUpdate)
        val subworkflowLabelsQuery = JsonQuery.compile("..|.subWorkflowMetadata? // empty|.labels|.short_id", Versions.JQ_1_5)
        subworkflowLabelsQuery.apply(scope, jsonNodeAfterUpdate, subworkflowOutputAfterUpdate)

        val actual = subworkflowNodesAfterUpdate.toList map { _.textValue() }
        val expected = subWorkflowIds map { _.shortString }
        actual shouldEqual expected

      case FailedMetadataResponse(_, reason) => fail(reason)
    }
  }
}

object CarbonitedMetadataThawingActorSpec {
  val workflowId = WorkflowId.fromString("2ce544a0-4c0d-4cc9-8a0b-b412bb1e5f82")

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
