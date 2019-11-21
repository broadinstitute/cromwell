package cromwell.engine.workflow.lifecycle.deletion

import akka.actor.ActorRef
import akka.testkit.{TestFSMRef, TestProbe}
import cromwell.core.path.Path
import cromwell.core.path.PathFactory.PathBuilders
import cromwell.core.{RootWorkflowId, TestKitSuite, WorkflowId}
import cromwell.engine.workflow.lifecycle.deletion.DeleteWorkflowFilesActor.{DeletingIntermediateFiles, DeletingIntermediateFilesData, FetchingAllOutputs, FetchingFinalOutputs, StartWorkflowFilesDeletion}
import cromwell.filesystems.gcs.{GcsPathBuilder, MockGcsPathBuilder}
import cromwell.services.metadata.MetadataService.{GetRootAndSubworkflowOutputs, RootAndSubworkflowOutputsLookupResponse, WorkflowOutputs, WorkflowOutputsResponse}
import cromwell.services.metadata._
import org.scalatest.concurrent.Eventually.eventually
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._

class DeleteWorkflowFilesActorSpec extends TestKitSuite("DeleteWorkflowFilesActorSpec") with FlatSpecLike with Matchers {

  val serviceRegistryActor = TestProbe()

  val mockPathBuilder: GcsPathBuilder = MockGcsPathBuilder.instance
  val mockPathBuilders = List(mockPathBuilder)

  it should "follow the expected golden-path lifecycle" in {
    val rootWorkflowId = RootWorkflowId(WorkflowId.randomId().id)
    val subworkflowId = WorkflowId.randomId()
    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, serviceRegistryActor.ref, mockPathBuilders))

    val allOutputsMetadataEvents = Vector(
      MetadataEvent(
        MetadataKey(subworkflowId, Option(MetadataJobKey("first_sub_workflow.first_task",None,1)), "outputs:first_task_output_1"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/output_file.txt", MetadataString))),
      MetadataEvent(
        MetadataKey(subworkflowId, Option(MetadataJobKey("first_sub_workflow.first_task",None,1)), "outputs:first_task_output_2"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/intermediate_file.txt", MetadataString))),
      MetadataEvent(
        MetadataKey(subworkflowId,None,"outputs:first_sub_workflow.second_output_file"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/intermediate_file.txt",MetadataString))),
      MetadataEvent(
        MetadataKey(subworkflowId,None,"outputs:first_sub_workflow.first_output_file"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/output_file.txt",MetadataString))),
      MetadataEvent(
        MetadataKey(rootWorkflowId,Option(MetadataJobKey("main_workflow.first_sub_workflow",None,1)),"outputs:second_output_file"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/intermediate_file.txt",MetadataString))),
      MetadataEvent(
        MetadataKey(rootWorkflowId,Option(MetadataJobKey("main_workflow.first_sub_workflow",None,1)),"outputs:first_output_file"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/output_file.txt",MetadataString))),
      MetadataEvent(
        MetadataKey(rootWorkflowId,None,"outputs:main_workflow.main_output"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/output_file.txt",MetadataString)))
    )

    val finalOutputsMetadataEvents = Vector(
      MetadataEvent(
        MetadataKey(rootWorkflowId,None,"outputs:main_workflow.main_output"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/output_file.txt",MetadataString)))
    )

    val gcsFilePath: Path = mockPathBuilder.build(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/intermediate_file.txt").get
    val intermediateFileList = Set(gcsFilePath)

    testDeleteWorkflowFilesActor ! StartWorkflowFilesDeletion

    serviceRegistryActor.expectMsgPF(10.seconds) {
      case m: GetRootAndSubworkflowOutputs => m.workflowId shouldBe rootWorkflowId
    }
    eventually {
      testDeleteWorkflowFilesActor.stateName shouldBe FetchingAllOutputs
    }

    serviceRegistryActor.send(testDeleteWorkflowFilesActor, RootAndSubworkflowOutputsLookupResponse(rootWorkflowId, allOutputsMetadataEvents))

    serviceRegistryActor.expectMsgPF(10.seconds) {
      case o: WorkflowOutputs =>
        o.workflowId shouldBe rootWorkflowId
        o.convertResponseToJson shouldBe false
    }
    eventually {
      testDeleteWorkflowFilesActor.stateName shouldBe FetchingFinalOutputs
    }

    serviceRegistryActor.send(testDeleteWorkflowFilesActor, WorkflowOutputsResponse(rootWorkflowId, finalOutputsMetadataEvents))

    eventually {
      testDeleteWorkflowFilesActor.stateName shouldBe DeletingIntermediateFiles
      testDeleteWorkflowFilesActor.stateData shouldBe DeletingIntermediateFilesData(intermediateFileList)
    }

    //TODO: add tests as needed for deleting files step as part of WA-41 (https://broadworkbench.atlassian.net/browse/WA-41)
  }


  it should "remove any non-file intermediate outputs " in {
    val rootWorkflowId = RootWorkflowId(WorkflowId.randomId().id)
    val subworkflowId = WorkflowId.randomId()
    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, serviceRegistryActor.ref, mockPathBuilders))

    val allOutputsMetadataEvents = Vector(
      MetadataEvent(
        MetadataKey(subworkflowId, Option(MetadataJobKey("first_sub_workflow.first_task",None,1)), "outputs:first_task_output_1"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/output_file.txt", MetadataString))),
      MetadataEvent(
        MetadataKey(subworkflowId, Option(MetadataJobKey("first_sub_workflow.first_task",None,1)), "outputs:first_task_output_2"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/intermediate_file.txt", MetadataString))),
      MetadataEvent(
        MetadataKey(subworkflowId,None,"outputs:first_sub_workflow.second_output_file"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/intermediate_file.txt",MetadataString))),
      MetadataEvent(
        MetadataKey(subworkflowId,None,"outputs:first_sub_workflow.first_output_file"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/output_file.txt",MetadataString))),
      MetadataEvent(
        MetadataKey(rootWorkflowId,Option(MetadataJobKey("main_workflow.first_sub_workflow",None,1)),"outputs:second_output_file"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/intermediate_file.txt",MetadataString))),
      MetadataEvent(
        MetadataKey(rootWorkflowId,Option(MetadataJobKey("main_workflow.first_sub_workflow",None,1)),"outputs:first_output_file"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/output_file.txt",MetadataString))),
      MetadataEvent(
        MetadataKey(rootWorkflowId,None,"outputs:main_workflow.main_output"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/output_file.txt",MetadataString))),
      MetadataEvent(
        MetadataKey(rootWorkflowId,Option(MetadataJobKey("main_workflow.first_sub_workflow",None,1)),"outputs:string_output"),
        Option(MetadataValue("Hello World!",MetadataString))),
      MetadataEvent(
        MetadataKey(rootWorkflowId,None,"outputs:main_workflow.main_string_output"),
        Option(MetadataValue("Hello World!",MetadataString))),
      MetadataEvent(
        MetadataKey(subworkflowId,None,"outputs:first_sub_workflow.string_output"),
        Option(MetadataValue("Hello World!",MetadataString)))
    )

    val finalOutputsMetadataEvents = Vector(
      MetadataEvent(
        MetadataKey(rootWorkflowId,None,"outputs:main_workflow.main_output"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/output_file.txt",MetadataString))),
      MetadataEvent(
        MetadataKey(rootWorkflowId,None,"outputs:main_workflow.main_string_output"),
        Option(MetadataValue("Hello World!",MetadataString)))
    )

    val gcsFilePath: Path = mockPathBuilder.build(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/intermediate_file.txt").get
    val expectedIntermediateFileList = Set(gcsFilePath)

    val actualIntermediateFileList = testDeleteWorkflowFilesActor.underlyingActor.gatherIntermediateOutputFiles(allOutputsMetadataEvents, finalOutputsMetadataEvents)

    actualIntermediateFileList shouldBe expectedIntermediateFileList
  }


  it should "terminate if root workflow has no outputs" in {
    val testProbe = TestProbe()
    val rootWorkflowId = RootWorkflowId(WorkflowId.randomId().id)
    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, serviceRegistryActor.ref, mockPathBuilders))

    testProbe watch testDeleteWorkflowFilesActor

    testDeleteWorkflowFilesActor ! StartWorkflowFilesDeletion

    serviceRegistryActor.expectMsgPF(10.seconds) {
      case m: GetRootAndSubworkflowOutputs => m.workflowId shouldBe rootWorkflowId
    }
    eventually {
      testDeleteWorkflowFilesActor.stateName shouldBe FetchingAllOutputs
    }

    serviceRegistryActor.send(testDeleteWorkflowFilesActor, RootAndSubworkflowOutputsLookupResponse(rootWorkflowId, Vector.empty[MetadataEvent]))

    testProbe.expectTerminated(testDeleteWorkflowFilesActor, 10.seconds)
  }


  it should "terminate if root workflow has no intermediate outputs to delete" in {
    val testProbe = TestProbe()
    val rootWorkflowId = RootWorkflowId(WorkflowId.randomId().id)
    val subworkflowId = WorkflowId.randomId()
    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, serviceRegistryActor.ref, mockPathBuilders))

    val allOutputsMetadataEvents = Vector(
      MetadataEvent(
        MetadataKey(subworkflowId, Option(MetadataJobKey("first_sub_workflow.first_task",None,1)), "outputs:first_task_output_1"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/output_file.txt", MetadataString))),
      MetadataEvent(
        MetadataKey(subworkflowId,None,"outputs:first_sub_workflow.first_output_file"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/output_file.txt",MetadataString))),
      MetadataEvent(
        MetadataKey(rootWorkflowId,Option(MetadataJobKey("main_workflow.first_sub_workflow",None,1)),"outputs:first_output_file"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/output_file.txt",MetadataString))),
      MetadataEvent(
        MetadataKey(rootWorkflowId,None,"outputs:main_workflow.main_output"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/output_file.txt",MetadataString)))
    )

    val finalOutputsMetadataEvents = Vector(
      MetadataEvent(
        MetadataKey(rootWorkflowId,None,"outputs:main_workflow.main_output"),
        Option(MetadataValue(s"gs://my-bucket/main_workflow/${rootWorkflowId.toString}/subworkflow-path/${subworkflowId.toString}/call-first_task/output_file.txt",MetadataString)))
    )

    testProbe watch testDeleteWorkflowFilesActor

    testDeleteWorkflowFilesActor ! StartWorkflowFilesDeletion

    serviceRegistryActor.expectMsgPF(10.seconds) {
      case m: GetRootAndSubworkflowOutputs => m.workflowId shouldBe rootWorkflowId
    }
    eventually {
      testDeleteWorkflowFilesActor.stateName shouldBe FetchingAllOutputs
    }

    serviceRegistryActor.send(testDeleteWorkflowFilesActor, RootAndSubworkflowOutputsLookupResponse(rootWorkflowId, allOutputsMetadataEvents))

    serviceRegistryActor.expectMsgPF(10.seconds) {
      case o: WorkflowOutputs =>
        o.workflowId shouldBe rootWorkflowId
        o.convertResponseToJson shouldBe false
    }
    eventually {
      testDeleteWorkflowFilesActor.stateName shouldBe FetchingFinalOutputs
    }

    serviceRegistryActor.send(testDeleteWorkflowFilesActor, WorkflowOutputsResponse(rootWorkflowId, finalOutputsMetadataEvents))

    testProbe.expectTerminated(testDeleteWorkflowFilesActor, 10.seconds)
  }
}


class MockDeleteWorkflowFilesActor(rootWorkflowId: RootWorkflowId, serviceRegistryActor: ActorRef, pathBuilders: PathBuilders) extends
  DeleteWorkflowFilesActor(rootWorkflowId, serviceRegistryActor, pathBuilders) { }
