package cromwell.engine.workflow.lifecycle.deletion

import java.io.FileNotFoundException

import akka.actor.ActorRef
import akka.testkit.{TestFSMRef, TestProbe}
import common.util.Backoff
import cromwell.core._
import cromwell.core.actor.StreamIntegration.BackPressure
import cromwell.core.io.{IoCommandBuilder, IoDeleteCommand, IoFailure, IoSuccess, PartialIoCommandBuilder}
import cromwell.core.path.Path
import cromwell.core.path.PathFactory.PathBuilders
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.engine.io.IoAttempts.EnhancedCromwellIoException
import cromwell.engine.workflow.lifecycle.deletion.DeleteWorkflowFilesActor.{StartWorkflowFilesDeletion, WaitingForIoResponses}
import cromwell.filesystems.gcs.batch.{GcsBatchCommandBuilder, GcsBatchDeleteCommand}
import cromwell.filesystems.gcs.{GcsPath, GcsPathBuilder, MockGcsPathBuilder}
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.impl.FileDeletionStatus
import cromwell.services.metadata.impl.FileDeletionStatus.{Failed, InProgress, Succeeded}
import org.scalatest.BeforeAndAfter
import org.scalatest.concurrent.Eventually.eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import wom.graph.GraphNodePort.{ExpressionBasedOutputPort, GraphNodeOutputPort}
import wom.graph.{FullyQualifiedName, LocalName, WomIdentifier}
import wom.types.{WomMaybeEmptyArrayType, WomSingleFileType, WomStringType}
import wom.values.{WomArray, WomSingleFile, WomString}

import scala.concurrent.duration._
import scala.util.control.NoStackTrace
import scala.util.{Failure, Try}

class DeleteWorkflowFilesActorSpec extends TestKitSuite
  with AnyFlatSpecLike
  with Matchers
  with BeforeAndAfter {

  val mockPathBuilder: GcsPathBuilder = MockGcsPathBuilder.instance
  val mockPathBuilders = List(mockPathBuilder)
  private val serviceRegistryActor = TestProbe("serviceRegistryActor")
  private val ioActor = TestProbe("ioActor")

  val emptyWorkflowIdSet = Set.empty[WorkflowId]

  var testProbe: TestProbe = _
  var rootWorkflowId: RootWorkflowId = _
  var subworkflowId: WorkflowId = _
  var allOutputs: CallOutputs = _
  var finalOutputs: CallOutputs = _
  var testDeleteWorkflowFilesActor: TestFSMRef[DeleteWorkflowFilesActor.DeleteWorkflowFilesActorState, DeleteWorkflowFilesActor.DeleteWorkflowFilesActorStateData, MockDeleteWorkflowFilesActor] = _
  var gcsFilePath: GcsPath = _

  before {
    rootWorkflowId = RootWorkflowId(WorkflowId.randomId().id)
    subworkflowId = WorkflowId.randomId()
    testProbe = TestProbe(s"test-probe-$rootWorkflowId")

    allOutputs = CallOutputs(Map(
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
    finalOutputs = CallOutputs(Map(
      GraphNodeOutputPort(WomIdentifier(LocalName("main_output"),FullyQualifiedName("main_workflow.main_output")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt")
    ))

    testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, emptyWorkflowIdSet, finalOutputs, allOutputs, mockPathBuilders, serviceRegistryActor.ref, ioActor.ref))
    testProbe watch testDeleteWorkflowFilesActor

    gcsFilePath = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt").get
  }

  it should "follow the expected golden-path lifecycle" in {

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

    testDeleteWorkflowFilesActor !
      IoSuccess(GcsBatchDeleteCommand.forPath(gcsFilePath, swallowIOExceptions = false).get, ())

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

  it should "be resilient to backpressure from the IoActor in an otherwise golden-path lifecycle" in {

    testDeleteWorkflowFilesActor ! StartWorkflowFilesDeletion

    serviceRegistryActor.expectMsgPF(10.seconds) {
      case m: PutMetadataAction =>
        val event = m.events.head
        m.events.size shouldBe 1
        event.key.workflowId shouldBe rootWorkflowId
        event.key.key shouldBe WorkflowMetadataKeys.FileDeletionStatus
        event.value.get.value shouldBe FileDeletionStatus.toDatabaseValue(InProgress)
    }

    val expectedDeleteCommand = GcsBatchDeleteCommand.forPath(gcsFilePath, swallowIOExceptions = false).get

    ioActor.expectMsgPF(10.seconds) {
      case `expectedDeleteCommand` => // woohoo!
    }

    eventually {
      testDeleteWorkflowFilesActor.stateName shouldBe WaitingForIoResponses
    }

    // Simulate a few backpressure events before we get the success:
    0 until 10 foreach { _ =>
      ioActor.send(testDeleteWorkflowFilesActor, BackPressure(expectedDeleteCommand))
      ioActor.expectMsgPF(10.seconds) {
        case cmd: IoDeleteCommand => cmd.file shouldBe gcsFilePath
        case other => fail(s"Bad message: $other")
      }
    }

    ioActor.send(testDeleteWorkflowFilesActor, IoSuccess(expectedDeleteCommand, ()))

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

    testDeleteWorkflowFilesActor !
      IoFailure(
        command = GcsBatchDeleteCommand.forPath(gcsFilePath, swallowIOExceptions = false).get,
        failure = new Exception(s"Something is fishy!"),
      )

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
    testDeleteWorkflowFilesActor !
      IoFailure(
        command = GcsBatchDeleteCommand.forPath(gcsFilePath, swallowIOExceptions = false).get,
        failure = fileNotFoundException,
      )

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

    val expectedIntermediateFiles = Set(gcsFilePath)

    val actualIntermediateFiles = testDeleteWorkflowFilesActor.underlyingActor.gatherIntermediateOutputFiles(
      allOutputs.outputs.values.toSet,
      finalOutputs.outputs.values.toSet
    )

    actualIntermediateFiles shouldBe expectedIntermediateFiles
  }


  it should "delete all file outputs if there are no final outputs" in {

    finalOutputs = CallOutputs.empty

    val gcsFilePath1: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt").get
    val gcsFilePath2: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt").get
    val expectedIntermediateFiles = Set(gcsFilePath1, gcsFilePath2)

    testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, emptyWorkflowIdSet, finalOutputs, allOutputs, mockPathBuilders, serviceRegistryActor.ref, ioActor.ref))

    val actualIntermediateFiles = testDeleteWorkflowFilesActor.underlyingActor.gatherIntermediateOutputFiles(
      allOutputs.outputs.values.toSet,
      finalOutputs.outputs.values.toSet
    )

    actualIntermediateFiles shouldBe expectedIntermediateFiles
  }

  it should "send failure when delete command creation is unsuccessful for a file" in {

    val partialIoCommandBuilder = new PartialIoCommandBuilder {
      override def deleteCommand: PartialFunction[(Path, Boolean), Try[IoDeleteCommand]] = {
        case _ => Failure(new Exception("everything's fine, I am an expected delete fail") with NoStackTrace)
      }
    }
    val ioCommandBuilder = new IoCommandBuilder(List(partialIoCommandBuilder))

    testDeleteWorkflowFilesActor =
      TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId,
        rootAndSubworkflowIds = emptyWorkflowIdSet,
        workflowFinalOutputs = finalOutputs,
        workflowAllOutputs = allOutputs,
        pathBuilders = mockPathBuilders,
        serviceRegistryActor = serviceRegistryActor.ref,
        ioActor = ioActor.ref,
        gcsCommandBuilder = ioCommandBuilder,
      ))
    testProbe.watch(testDeleteWorkflowFilesActor)

    testDeleteWorkflowFilesActor ! StartWorkflowFilesDeletion

    serviceRegistryActor.expectMsgPF(10.seconds) {
      case m: PutMetadataAction =>
        val event = m.events.head
        m.events.size shouldBe 1
        event.key.workflowId shouldBe rootWorkflowId
        event.key.key shouldBe WorkflowMetadataKeys.FileDeletionStatus
        event.value.get.value shouldBe FileDeletionStatus.toDatabaseValue(InProgress)
    }

    testProbe.expectTerminated(testDeleteWorkflowFilesActor, 10.seconds)
    testProbe.unwatch(testDeleteWorkflowFilesActor)
  }

  it should "terminate if root workflow has no intermediate outputs to delete" in {

    finalOutputs = CallOutputs(Map(
      GraphNodeOutputPort(WomIdentifier(LocalName("main_output_1"),FullyQualifiedName("main_workflow.main_output_1")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt") ,
      GraphNodeOutputPort(WomIdentifier(LocalName("main_output_2"),FullyQualifiedName("main_workflow.main_output_2")), WomSingleFileType, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file2.txt")
    ))

    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, emptyWorkflowIdSet, finalOutputs, allOutputs, mockPathBuilders, serviceRegistryActor.ref, ioActor.ref))

    testProbe watch testDeleteWorkflowFilesActor

    testDeleteWorkflowFilesActor ! StartWorkflowFilesDeletion

    testProbe.expectTerminated(testDeleteWorkflowFilesActor, 10.seconds)
  }


  it should "remove values that are file names in form of string" in {

    allOutputs = allOutputs.copy(outputs = allOutputs.outputs ++ Map(
      GraphNodeOutputPort(WomIdentifier(LocalName("file_with_file_path_string_output"),FullyQualifiedName("first_sub_workflow.file_with_file_path_string_output")), WomStringType, null)
        -> WomString(s"gs://my_bucket/non_existent_file.txt"),
      ExpressionBasedOutputPort(WomIdentifier(LocalName("first_task.first_task_file_with_file_path_output"),FullyQualifiedName("first_sub_workflow.first_task.first_task_file_with_file_path_output")), WomSingleFileType, null, null)
        -> WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/file_with_file_path.txt")
    ))

    val gcsFilePath1: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt").get
    val gcsFilePath2: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/file_with_file_path.txt").get
    val expectedIntermediateFiles = Set(gcsFilePath1, gcsFilePath2)

    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, emptyWorkflowIdSet, finalOutputs, allOutputs, mockPathBuilders, serviceRegistryActor.ref, ioActor.ref))

    val actualIntermediateFiles = testDeleteWorkflowFilesActor.underlyingActor.gatherIntermediateOutputFiles(
      allOutputs.outputs.values.toSet,
      finalOutputs.outputs.values.toSet
    )

    actualIntermediateFiles shouldBe expectedIntermediateFiles
  }


  it should "identify and gather glob files" in {

    allOutputs = allOutputs.copy(outputs = allOutputs.outputs ++ Map(
      GraphNodeOutputPort(WomIdentifier(LocalName("glob_output"),FullyQualifiedName("first_sub_workflow.glob_output")), WomMaybeEmptyArrayType(WomSingleFileType), null)
        -> WomArray(Seq(WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/glob-random_id/intermediate_file1.txt"),
        WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/glob-random_id/intermediate_file2.txt"))),
      ExpressionBasedOutputPort(WomIdentifier(LocalName("first_task.first_task_glob"),FullyQualifiedName("first_sub_workflow.first_task.first_task_glob")), WomMaybeEmptyArrayType(WomSingleFileType), null, null)
        -> WomArray(Seq(WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/glob-random_id/intermediate_file1.txt"),
        WomSingleFile(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/glob-random_id/intermediate_file2.txt")))
    ))

    val gcsFilePath1: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/intermediate_file1.txt").get
    val gcsGlobFilePath1: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/glob-random_id/intermediate_file1.txt").get
    val gcsGlobFilePath2: Path = mockPathBuilder.build(s"gs://my_bucket/main_workflow/$rootWorkflowId/call-first_sub_workflow/firstSubWf.first_sub_workflow/$subworkflowId/call-first_task/glob-random_id/intermediate_file2.txt").get

    val expectedIntermediateFiles = Set(gcsFilePath1, gcsGlobFilePath1, gcsGlobFilePath2)

    val testDeleteWorkflowFilesActor = TestFSMRef(new MockDeleteWorkflowFilesActor(rootWorkflowId, emptyWorkflowIdSet, finalOutputs, allOutputs, mockPathBuilders, serviceRegistryActor.ref, ioActor.ref))

    val actualIntermediateFiles = testDeleteWorkflowFilesActor.underlyingActor.gatherIntermediateOutputFiles(
      allOutputs.outputs.values.toSet,
      finalOutputs.outputs.values.toSet
    )

    actualIntermediateFiles shouldBe expectedIntermediateFiles
  }
}


class MockDeleteWorkflowFilesActor(rootWorkflowId: RootWorkflowId,
                                   rootAndSubworkflowIds: Set[WorkflowId],
                                   workflowFinalOutputs: CallOutputs,
                                   workflowAllOutputs: CallOutputs,
                                   pathBuilders: PathBuilders,
                                   serviceRegistryActor: ActorRef,
                                   ioActor: ActorRef,
                                   gcsCommandBuilder: IoCommandBuilder = GcsBatchCommandBuilder,
                                  ) extends
  DeleteWorkflowFilesActor(
    rootWorkflowId,
    rootAndSubworkflowIds,
    workflowFinalOutputs.outputs.values.toSet,
    workflowAllOutputs.outputs.values.toSet,
    pathBuilders,
    serviceRegistryActor,
    ioActor,
    gcsCommandBuilder,
  ) {
  // Override the IO actor backoff for the benefit of the backpressure tests:
  override def initialBackoff(): Backoff = SimpleExponentialBackoff(100.millis, 1.second, 1.2D)
}
