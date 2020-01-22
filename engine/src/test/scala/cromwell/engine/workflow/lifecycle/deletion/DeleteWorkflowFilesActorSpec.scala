package cromwell.engine.workflow.lifecycle.deletion

import java.io.FileNotFoundException

import akka.actor.ActorRef
import akka.testkit.{TestFSMRef, TestProbe}
import cromwell.core._
import cromwell.core.io.{IoDeleteCommand, IoFailure, IoSuccess}
import cromwell.core.path.Path
import cromwell.core.path.PathFactory.PathBuilders
import cromwell.engine.io.IoAttempts.EnhancedCromwellIoException
import cromwell.engine.workflow.lifecycle.deletion.DeleteWorkflowFilesActor.{StartWorkflowFilesDeletion, WaitingForIoResponses}
import cromwell.filesystems.gcs.batch.GcsBatchDeleteCommand
import cromwell.filesystems.gcs.{GcsPath, GcsPathBuilder, MockGcsPathBuilder}
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.impl.FileDeletionStatus
import cromwell.services.metadata.impl.FileDeletionStatus.{Failed, InProgress, Succeeded}
import org.scalatest.concurrent.Eventually.eventually
import org.scalatest.{FlatSpecLike, Matchers}
import wom.graph.GraphNodePort.{ExpressionBasedOutputPort, GraphNodeOutputPort}
import wom.graph.{FullyQualifiedName, LocalName, WomIdentifier}
import wom.types.{WomMaybeEmptyArrayType, WomSingleFileType, WomStringType}
import wom.values.{WomArray, WomSingleFile, WomString}

import scala.concurrent.duration._

class DeleteWorkflowFilesActorSpec extends TestKitSuite("DeleteWorkflowFilesActorSpec") with FlatSpecLike with Matchers {

  val mockPathBuilder: GcsPathBuilder = MockGcsPathBuilder.instance
  val mockPathBuilders = List(mockPathBuilder)
  val serviceRegistryActor = TestProbe()
  val ioActor = TestProbe()

  it should "follow the expected golden-path lifecycle" in {
    val testProbe = TestProbe()
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

    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, finalOutputs, allOutputs, mockPathBuilders, serviceRegistryActor.ref, ioActor.ref))
    testProbe watch testDeleteWorkflowFilesActor

    val gcsFilePath: GcsPath = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt").get

    testDeleteWorkflowFilesActor ! StartWorkflowFilesDeletion

    serviceRegistryActor.expectMsgPF(10.seconds) {
      case m: PutMetadataAction =>
        val event = m.events.head
        m.events.size shouldBe 1
        event.key.workflowId shouldBe rootWorkflowId
        event.key.key shouldBe WorkflowMetadataKeys.FileDeletionStatus
        event.value.get.value shouldBe FileDeletionStatus.toDatabaseValue(InProgress)
    }

    ioActor.expectMsgPF(10.seconds) {
      case cmd: IoDeleteCommand => cmd.file shouldBe gcsFilePath
    }

    eventually {
      testDeleteWorkflowFilesActor.stateName shouldBe WaitingForIoResponses
    }

    testDeleteWorkflowFilesActor ! IoSuccess(GcsBatchDeleteCommand(gcsFilePath, swallowIOExceptions = false), ())

    serviceRegistryActor.expectMsgPF(10.seconds) {
      case m: PutMetadataAction =>
        val event = m.events.head
        m.events.size shouldBe 1
        event.key.workflowId shouldBe rootWorkflowId
        event.key.key shouldBe WorkflowMetadataKeys.FileDeletionStatus
        event.value.get.value shouldBe FileDeletionStatus.toDatabaseValue(Succeeded)
    }

    testProbe.expectTerminated(testDeleteWorkflowFilesActor, 10.seconds)
  }


  it should "send failure when delete is unsuccessful for a file" in {
    val testProbe = TestProbe()
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

    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, finalOutputs, allOutputs, mockPathBuilders, serviceRegistryActor.ref, ioActor.ref))
    testProbe watch testDeleteWorkflowFilesActor

    val gcsFilePath: GcsPath = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt").get

    testDeleteWorkflowFilesActor ! StartWorkflowFilesDeletion

    serviceRegistryActor.expectMsgPF(10.seconds) {
      case m: PutMetadataAction =>
        val event = m.events.head
        m.events.size shouldBe 1
        event.key.workflowId shouldBe rootWorkflowId
        event.key.key shouldBe WorkflowMetadataKeys.FileDeletionStatus
        event.value.get.value shouldBe FileDeletionStatus.toDatabaseValue(InProgress)
    }

    ioActor.expectMsgPF(10.seconds) {
      case cmd: IoDeleteCommand => cmd.file shouldBe gcsFilePath
    }

    eventually {
      testDeleteWorkflowFilesActor.stateName shouldBe WaitingForIoResponses
    }

    testDeleteWorkflowFilesActor ! IoFailure(GcsBatchDeleteCommand(gcsFilePath, swallowIOExceptions = false), new Exception(s"Something is fishy!"))

    serviceRegistryActor.expectMsgPF(10.seconds) {
      case m: PutMetadataAction =>
        val event = m.events.head
        m.events.size shouldBe 1
        event.key.workflowId shouldBe rootWorkflowId
        event.key.key shouldBe WorkflowMetadataKeys.FileDeletionStatus
        event.value.get.value shouldBe FileDeletionStatus.toDatabaseValue(Failed)
    }

    testProbe.expectTerminated(testDeleteWorkflowFilesActor, 10.seconds)
  }


  it should "send success when the failure is FileNotFound" in {
    val testProbe = TestProbe()
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

    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, finalOutputs, allOutputs, mockPathBuilders, serviceRegistryActor.ref, ioActor.ref))
    testProbe watch testDeleteWorkflowFilesActor

    val gcsFilePath: GcsPath = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt").get

    testDeleteWorkflowFilesActor ! StartWorkflowFilesDeletion

    serviceRegistryActor.expectMsgPF(10.seconds) {
      case m: PutMetadataAction =>
        val event = m.events.head
        m.events.size shouldBe 1
        event.key.workflowId shouldBe rootWorkflowId
        event.key.key shouldBe WorkflowMetadataKeys.FileDeletionStatus
        event.value.get.value shouldBe FileDeletionStatus.toDatabaseValue(InProgress)
    }

    ioActor.expectMsgPF(10.seconds) {
      case cmd: IoDeleteCommand => cmd.file shouldBe gcsFilePath
    }

    eventually {
      testDeleteWorkflowFilesActor.stateName shouldBe WaitingForIoResponses
    }

    val fileNotFoundException = EnhancedCromwellIoException(s"File not found", new FileNotFoundException(gcsFilePath.pathAsString))
    testDeleteWorkflowFilesActor ! IoFailure(GcsBatchDeleteCommand(gcsFilePath, swallowIOExceptions = false), fileNotFoundException)

    serviceRegistryActor.expectMsgPF(10.seconds) {
      case m: PutMetadataAction =>
        val event = m.events.head
        m.events.size shouldBe 1
        event.key.workflowId shouldBe rootWorkflowId
        event.key.key shouldBe WorkflowMetadataKeys.FileDeletionStatus
        event.value.get.value shouldBe FileDeletionStatus.toDatabaseValue(Succeeded)
    }

    testProbe.expectTerminated(testDeleteWorkflowFilesActor, 10.seconds)
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

    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, finalOutputs, allOutputs, mockPathBuilders, serviceRegistryActor.ref, ioActor.ref))

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

    val gcsFilePath1: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt").get
    val gcsFilePath2: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt").get
    val expectedIntermediateFileList = Set(gcsFilePath1, gcsFilePath2)

    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, finalOutputs, allOutputs, mockPathBuilders, serviceRegistryActor.ref, ioActor.ref))

    val actualIntermediateFileList = testDeleteWorkflowFilesActor.underlyingActor.gatherIntermediateOutputFiles(allOutputs.outputs, finalOutputs.outputs)

    actualIntermediateFileList shouldBe expectedIntermediateFileList
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

    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, finalOutputs, allOutputs, mockPathBuilders, serviceRegistryActor.ref, ioActor.ref))

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

    val gcsFilePath1: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt").get
    val gcsFilePath2: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/file_with_file_path.txt").get
    val expectedIntermediateFileList = Set(gcsFilePath1, gcsFilePath2)

    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, finalOutputs, allOutputs, mockPathBuilders, serviceRegistryActor.ref, ioActor.ref))

    val actualIntermediateFileList = testDeleteWorkflowFilesActor.underlyingActor.gatherIntermediateOutputFiles(allOutputs.outputs, finalOutputs.outputs)

    actualIntermediateFileList shouldBe expectedIntermediateFileList
  }


  it should "identify and gather glob files" in {
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
      GraphNodeOutputPort(WomIdentifier(LocalName("glob_output"),FullyQualifiedName("first_sub_workflow.glob_output")), WomMaybeEmptyArrayType(WomSingleFileType), null)
        -> WomArray(Seq(WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/glob-random_id/intermediate_file1.txt"),
        WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/glob-random_id/intermediate_file2.txt"))),
      ExpressionBasedOutputPort(WomIdentifier(LocalName("first_task.first_task_glob"),FullyQualifiedName("first_sub_workflow.first_task.first_task_glob")), WomMaybeEmptyArrayType(WomSingleFileType), null, null)
        -> WomArray(Seq(WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/glob-random_id/intermediate_file1.txt"),
        WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/glob-random_id/intermediate_file2.txt")))
    ))

    val finalOutputs = CallOutputs(Map(
      GraphNodeOutputPort(WomIdentifier(LocalName("main_output"),FullyQualifiedName("main_workflow.main_output")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt")
    ))

    val gcsFilePath1: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt").get
    val gcsGlobFilePath1: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/glob-random_id/intermediate_file1.txt").get
    val gcsGlobFilePath2: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/glob-random_id/intermediate_file2.txt").get

    val expectedIntermediateFileList = Set(gcsFilePath1, gcsGlobFilePath1, gcsGlobFilePath2)

    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, finalOutputs, allOutputs, mockPathBuilders, serviceRegistryActor.ref, ioActor.ref))

    val actualIntermediateFileList = testDeleteWorkflowFilesActor.underlyingActor.gatherIntermediateOutputFiles(allOutputs.outputs, finalOutputs.outputs)

    actualIntermediateFileList shouldBe expectedIntermediateFileList
  }
}


class MockDeleteWorkflowFilesActor(rootWorkflowId: RootWorkflowId, workflowFinalOutputs: CallOutputs,
                                   workflowAllOutputs: CallOutputs, pathBuilders: PathBuilders,
                                   serviceRegistryActor: ActorRef, ioActor: ActorRef) extends
  DeleteWorkflowFilesActor(rootWorkflowId, workflowFinalOutputs, workflowAllOutputs, pathBuilders, serviceRegistryActor, ioActor) { }
