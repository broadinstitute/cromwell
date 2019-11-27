package cromwell.engine.workflow.lifecycle.deletion

import akka.testkit.{TestFSMRef, TestProbe}
import cromwell.core.path.Path
import cromwell.core.path.PathFactory.PathBuilders
import cromwell.core.{CallOutputs, RootWorkflowId, TestKitSuite, WorkflowId}
import cromwell.engine.workflow.lifecycle.deletion.DeleteWorkflowFilesActor.{DeletingIntermediateFiles, DeletingIntermediateFilesData, StartWorkflowFilesDeletion}
import cromwell.filesystems.gcs.{GcsPathBuilder, MockGcsPathBuilder}
import org.scalatest.concurrent.Eventually.eventually
import org.scalatest.{FlatSpecLike, Matchers}
import wom.graph.GraphNodePort.{ExpressionBasedOutputPort, GraphNodeOutputPort}
import wom.graph.{FullyQualifiedName, LocalName, WomIdentifier}
import wom.types.{WomSingleFileType, WomStringType}
import wom.values.{WomSingleFile, WomString}

import scala.concurrent.duration._

class DeleteWorkflowFilesActorSpec extends TestKitSuite("DeleteWorkflowFilesActorSpec") with FlatSpecLike with Matchers {

  val mockPathBuilder: GcsPathBuilder = MockGcsPathBuilder.instance
  val mockPathBuilders = List(mockPathBuilder)

  it should "follow the expected golden-path lifecycle" in {
    val rootWorkflowId = RootWorkflowId(WorkflowId.randomId().id)
    val subworkflowId = WorkflowId.randomId()

    val allOutputs = CallOutputs(Map(
      GraphNodeOutputPort(WomIdentifier(LocalName("main_output"),FullyQualifiedName("main_workflow.main_output")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt"),
      ExpressionBasedOutputPort(WomIdentifier(LocalName("first_task.first_task_output_2"),FullyQualifiedName("first_sub_workflow.first_task.first_task_output_2")), WomSingleFileType, null, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt"),
      GraphNodeOutputPort(WomIdentifier(LocalName("first_output_file"),FullyQualifiedName("first_sub_workflow.first_output_file")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt"),
      ExpressionBasedOutputPort(WomIdentifier(LocalName("first_task.first_task_output_1"),FullyQualifiedName("first_sub_workflow.first_task.first_task_output_1")), WomSingleFileType, null, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt"),
      GraphNodeOutputPort(WomIdentifier(LocalName("second_output_file"),FullyQualifiedName("first_sub_workflow.second_output_file")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt")
    ))

    val finalOutputs = CallOutputs(Map(
      GraphNodeOutputPort(WomIdentifier(LocalName("main_output"),FullyQualifiedName("main_workflow.main_output")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt")
    ))

    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, finalOutputs, allOutputs, mockPathBuilders))

    val gcsFilePath: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt").get
    val intermediateFileList = Set(gcsFilePath)

    testDeleteWorkflowFilesActor ! StartWorkflowFilesDeletion

    eventually {
      testDeleteWorkflowFilesActor.stateName shouldBe DeletingIntermediateFiles
      testDeleteWorkflowFilesActor.stateData shouldBe DeletingIntermediateFilesData(intermediateFileList)
    }

    //TODO: add tests as needed for deleting files step as part of WA-41 (https://broadworkbench.atlassian.net/browse/WA-41)
  }


  it should "remove any non-file intermediate outputs" in {
    val rootWorkflowId = RootWorkflowId(WorkflowId.randomId().id)
    val subworkflowId = WorkflowId.randomId()

    val allOutputs = CallOutputs(Map(
      GraphNodeOutputPort(WomIdentifier(LocalName("main_output"),FullyQualifiedName("main_workflow.main_output")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt"),
      ExpressionBasedOutputPort(WomIdentifier(LocalName("first_task.first_task_output_2"),FullyQualifiedName("first_sub_workflow.first_task.first_task_output_2")), WomSingleFileType, null, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt"),
      GraphNodeOutputPort(WomIdentifier(LocalName("first_output_file"),FullyQualifiedName("first_sub_workflow.first_output_file")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt"),
      ExpressionBasedOutputPort(WomIdentifier(LocalName("first_task.first_task_output_1"),FullyQualifiedName("first_sub_workflow.first_task.first_task_output_1")), WomSingleFileType, null, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt"),
      GraphNodeOutputPort(WomIdentifier(LocalName("second_output_file"),FullyQualifiedName("first_sub_workflow.second_output_file")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt") ,
      GraphNodeOutputPort(WomIdentifier(LocalName("main_string_output"),FullyQualifiedName("main_workflow.main_string_output")), WomStringType, null)
        -> WomString("Hello World!"),
      GraphNodeOutputPort(WomIdentifier(LocalName("string_output"),FullyQualifiedName("first_sub_workflow.string_output")), WomStringType, null)
        -> WomString("Hello World!")
    ))

    val finalOutputs = CallOutputs(Map(
      GraphNodeOutputPort(WomIdentifier(LocalName("main_output"),FullyQualifiedName("main_workflow.main_output")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt"),
      GraphNodeOutputPort(WomIdentifier(LocalName("main_string_output"),FullyQualifiedName("main_workflow.main_string_output")), WomStringType, null)
        -> WomString("Hello World!")
    ))

    val gcsFilePath: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt").get
    val expectedIntermediateFileList = Set(gcsFilePath)

    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, finalOutputs, allOutputs, mockPathBuilders))

    val actualIntermediateFileList = testDeleteWorkflowFilesActor.underlyingActor.gatherIntermediateOutputFiles(allOutputs.outputs, finalOutputs.outputs)

    actualIntermediateFileList shouldBe expectedIntermediateFileList
  }


  it should "delete all file outputs if there are no final outputs" in {
    val rootWorkflowId = RootWorkflowId(WorkflowId.randomId().id)
    val subworkflowId = WorkflowId.randomId()

    val allOutputs = CallOutputs(Map(
      GraphNodeOutputPort(WomIdentifier(LocalName("main_output"),FullyQualifiedName("main_workflow.main_output")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt"),
      ExpressionBasedOutputPort(WomIdentifier(LocalName("first_task.first_task_output_2"),FullyQualifiedName("first_sub_workflow.first_task.first_task_output_2")), WomSingleFileType, null, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt"),
      GraphNodeOutputPort(WomIdentifier(LocalName("first_output_file"),FullyQualifiedName("first_sub_workflow.first_output_file")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt"),
      ExpressionBasedOutputPort(WomIdentifier(LocalName("first_task.first_task_output_1"),FullyQualifiedName("first_sub_workflow.first_task.first_task_output_1")), WomSingleFileType, null, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt"),
      GraphNodeOutputPort(WomIdentifier(LocalName("second_output_file"),FullyQualifiedName("first_sub_workflow.second_output_file")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt")
    ))

    val finalOutputs = CallOutputs.empty

    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, finalOutputs, allOutputs, mockPathBuilders))

    val gcsFilePath1: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt").get
    val gcsFilePath2: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt").get
    val intermediateFileList = Set(gcsFilePath1, gcsFilePath2)

    testDeleteWorkflowFilesActor ! StartWorkflowFilesDeletion

    eventually {
      testDeleteWorkflowFilesActor.stateName shouldBe DeletingIntermediateFiles
      testDeleteWorkflowFilesActor.stateData shouldBe DeletingIntermediateFilesData(intermediateFileList)
    }
  }


  it should "terminate if root workflow has no intermediate outputs to delete" in {
    val testProbe = TestProbe()
    val rootWorkflowId = RootWorkflowId(WorkflowId.randomId().id)
    val subworkflowId = WorkflowId.randomId()

    val allOutputs = CallOutputs(Map(
      GraphNodeOutputPort(WomIdentifier(LocalName("main_output_1"),FullyQualifiedName("main_workflow.main_output_1")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt"),
      ExpressionBasedOutputPort(WomIdentifier(LocalName("first_task.first_task_output_2"),FullyQualifiedName("first_sub_workflow.first_task.first_task_output_2")), WomSingleFileType, null, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt"),
      GraphNodeOutputPort(WomIdentifier(LocalName("first_output_file"),FullyQualifiedName("first_sub_workflow.first_output_file")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt"),
      ExpressionBasedOutputPort(WomIdentifier(LocalName("first_task.first_task_output_1"),FullyQualifiedName("first_sub_workflow.first_task.first_task_output_1")), WomSingleFileType, null, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt"),
      GraphNodeOutputPort(WomIdentifier(LocalName("second_output_file"),FullyQualifiedName("first_sub_workflow.second_output_file")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt") ,
      GraphNodeOutputPort(WomIdentifier(LocalName("main_output_2"),FullyQualifiedName("main_workflow.main_output_2")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt")
    ))

    val finalOutputs = CallOutputs(Map(
      GraphNodeOutputPort(WomIdentifier(LocalName("main_output_1"),FullyQualifiedName("main_workflow.main_output_1")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt") ,
      GraphNodeOutputPort(WomIdentifier(LocalName("main_output_2"),FullyQualifiedName("main_workflow.main_output_2")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt")
    ))

    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, finalOutputs, allOutputs, mockPathBuilders))


    testProbe watch testDeleteWorkflowFilesActor

    testDeleteWorkflowFilesActor ! StartWorkflowFilesDeletion

    testProbe.expectTerminated(testDeleteWorkflowFilesActor, 10.seconds)
  }


  it should "remove values that are file names in form of string" in {
    val rootWorkflowId = RootWorkflowId(WorkflowId.randomId().id)
        val subworkflowId = WorkflowId.randomId()

    val allOutputs = CallOutputs(Map(
      GraphNodeOutputPort(WomIdentifier(LocalName("main_output"),FullyQualifiedName("main_workflow.main_output")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt"),
      ExpressionBasedOutputPort(WomIdentifier(LocalName("first_task.first_task_output_2"),FullyQualifiedName("first_sub_workflow.first_task.first_task_output_2")), WomSingleFileType, null, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt"),
      GraphNodeOutputPort(WomIdentifier(LocalName("first_output_file"),FullyQualifiedName("first_sub_workflow.first_output_file")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt"),
      ExpressionBasedOutputPort(WomIdentifier(LocalName("first_task.first_task_output_1"),FullyQualifiedName("first_sub_workflow.first_task.first_task_output_1")), WomSingleFileType, null, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt"),
      GraphNodeOutputPort(WomIdentifier(LocalName("second_output_file"),FullyQualifiedName("first_sub_workflow.second_output_file")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt"),
      GraphNodeOutputPort(WomIdentifier(LocalName("file_with_file_path_string_output"),FullyQualifiedName("first_sub_workflow.file_with_file_path_string_output")), WomStringType, null)
        -> WomString(s"gs://my_bucket/non_existent_file.txt"),
      ExpressionBasedOutputPort(WomIdentifier(LocalName("first_task.first_task_file_with_file_path_output"),FullyQualifiedName("first_sub_workflow.first_task.first_task_file_with_file_path_output")), WomSingleFileType, null, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/file_with_file_path.txt")
    ))

    val finalOutputs = CallOutputs(Map(
      GraphNodeOutputPort(WomIdentifier(LocalName("main_output"),FullyQualifiedName("main_workflow.main_output")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt")
    ))

    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, finalOutputs, allOutputs, mockPathBuilders))

    val gcsFilePath: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt").get
    val intermediateFileList = Set(gcsFilePath)

    testDeleteWorkflowFilesActor ! StartWorkflowFilesDeletion

    eventually {
      testDeleteWorkflowFilesActor.stateName shouldBe DeletingIntermediateFiles
      testDeleteWorkflowFilesActor.stateData shouldBe DeletingIntermediateFilesData(intermediateFileList)
    }
  }
}


class MockDeleteWorkflowFilesActor(rootWorkflowId: RootWorkflowId, workflowFinalOutputs: CallOutputs,
                                   workflowAllOutputs: CallOutputs, pathBuilders: PathBuilders) extends
  DeleteWorkflowFilesActor(rootWorkflowId, workflowFinalOutputs, workflowAllOutputs, pathBuilders) { }
