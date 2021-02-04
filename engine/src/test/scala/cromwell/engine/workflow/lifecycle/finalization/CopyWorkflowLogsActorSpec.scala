package cromwell.engine.workflow.lifecycle.finalization

import akka.testkit._
import cromwell.core.io.DefaultIoCommand.{DefaultIoCopyCommand, DefaultIoDeleteCommand}
import cromwell.core.io._
import cromwell.core.logging.WorkflowLogger
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.util.GracefulShutdownHelper
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import scala.concurrent.duration._
import scala.util.control.NoStackTrace
import scala.util.{Failure, Try}

class CopyWorkflowLogsActorSpec
  extends TestKitSuite with AnyFlatSpecLike with Matchers {

  behavior of "CopyWorkflowLogsActor"

  private val msgWait = 10.second.dilated
  private val serviceRegistryActor = TestProbe("testServiceRegistryActor")
  private val ioActor = TestProbe("testIoActor")
  private val deathWatch = TestProbe("deathWatch")
  private val tempDir = DefaultPathBuilder.createTempDirectory("tempDir.")

  override protected def beforeAll(): Unit = {
    super.beforeAll()
  }

  override protected def afterAll(): Unit = {
    super.afterAll()
    tempDir.delete(true)
  }

  it should "send a command to delete" in {
    val workflowId = WorkflowId.randomId()
    val destinationDir = tempDir.resolve("deleteAFile")
    val props =
      CopyWorkflowLogsActor.props(
        serviceRegistryActor = serviceRegistryActor.ref,
        ioActor = ioActor.ref,
        workflowLogConfigurationOption = WorkflowLogger.workflowLogConfiguration,
        copyCommandBuilder = DefaultIoCommandBuilder,
        deleteCommandBuilder = DefaultIoCommandBuilder,
      )
    val copyWorkflowLogsActor = system.actorOf(props, "testCopyWorkflowLogsActor")

    // The delete command uses singletons so we expect that directory to be used
    lazy val workflowLogPath = WorkflowLogger.workflowLogConfiguration.get.dir.resolve(s"workflow.$workflowId.log")
    copyWorkflowLogsActor ! CopyWorkflowLogsActor.Copy(workflowId, destinationDir)

    // We expect a copy command to be created and sent to the ioActor as a Tuple2
    val copyCommand =
      DefaultIoCopyCommand(workflowLogPath, destinationDir.resolve(workflowLogPath.name))
    ioActor.expectMsg(msgWait, (workflowId, copyCommand))

    // Tell the copyWorkflowLogsActor the copy failed as a Tuple2
    copyWorkflowLogsActor !
      ((
        workflowId,
        IoFailure(copyCommand, new Exception("everything's fine, I am an expected copy fail") with NoStackTrace),
      ))

    // There should now be a delete command sent to the ioActor
    val deleteCommand = DefaultIoDeleteCommand(workflowLogPath, swallowIOExceptions = true)
    ioActor.expectMsg(msgWait, deleteCommand)

    // Tell the copyWorkflowLogsActor the delete failed
    copyWorkflowLogsActor !
      IoFailure(deleteCommand, new Exception("everything's fine, I am an expected delete fail") with NoStackTrace)

    // Send a shutdown after the delete
    deathWatch.watch(copyWorkflowLogsActor)
    copyWorkflowLogsActor ! GracefulShutdownHelper.ShutdownCommand

    // Then the actor should shutdown
    deathWatch.expectTerminated(copyWorkflowLogsActor, msgWait)
  }

  it should "skip copy and still delete if unable to create a copy command" in {
    val workflowId = WorkflowId.randomId()
    val destinationPath = DefaultPathBuilder.createTempFile(s"test_file_$workflowId.", ".file", Option(tempDir))
    val partialIoCommandBuilder = new PartialIoCommandBuilder {
      override def copyCommand: PartialFunction[(Path, Path), Try[IoCopyCommand]] = {
        case _ => Failure(new Exception("everything's fine, I am an expected copy fail") with NoStackTrace)
      }
    }
    val ioCommandBuilder = new IoCommandBuilder(List(partialIoCommandBuilder))
    val props =
      CopyWorkflowLogsActor.props(
        serviceRegistryActor = serviceRegistryActor.ref,
        ioActor = ioActor.ref,
        workflowLogConfigurationOption = WorkflowLogger.workflowLogConfiguration,
        copyCommandBuilder = ioCommandBuilder,
      )
    val copyWorkflowLogsActor = system.actorOf(props, "testCopyWorkflowLogsActorFailCopy")

    // The delete command uses singletons so we expect that directory to be used
    lazy val workflowLogPath = WorkflowLogger.workflowLogConfiguration.get.dir.resolve(s"workflow.$workflowId.log")
    copyWorkflowLogsActor ! CopyWorkflowLogsActor.Copy(workflowId, destinationPath)

    // Because the copy failed we expect a delete command to be sent immediately to the ioActor
    val deleteCommand = DefaultIoDeleteCommand(workflowLogPath, swallowIOExceptions = true)
    ioActor.expectMsg(msgWait, deleteCommand)

    copyWorkflowLogsActor ! IoFailure(deleteCommand, new Exception("everything's fine, I am an expected delete fail") with NoStackTrace)

    // Send a shutdown after the delete
    deathWatch.watch(copyWorkflowLogsActor)
    copyWorkflowLogsActor ! GracefulShutdownHelper.ShutdownCommand

    // Then the actor should shutdown
    deathWatch.expectTerminated(copyWorkflowLogsActor, msgWait)
  }

  it should "skip copy and still delete if unable to create a copy command while shutting down" in {
    val workflowId = WorkflowId.randomId()
    val destinationPath = DefaultPathBuilder.createTempFile(s"test_file_$workflowId.", ".file", Option(tempDir))
    val partialIoCommandBuilder = new PartialIoCommandBuilder {
      override def copyCommand: PartialFunction[(Path, Path), Try[IoCopyCommand]] = {
        case _ => Failure(new Exception("everything's fine, I am an expected copy fail") with NoStackTrace)
      }
    }
    val ioCommandBuilder = new IoCommandBuilder(List(partialIoCommandBuilder))
    val props = CopyWorkflowLogsActor.props(
      serviceRegistryActor = serviceRegistryActor.ref,
      ioActor = ioActor.ref,
      workflowLogConfigurationOption = WorkflowLogger.workflowLogConfiguration,
      copyCommandBuilder = ioCommandBuilder,
    )
    val copyWorkflowLogsActor = system.actorOf(props, "testCopyWorkflowLogsActorFailCopyShutdown")

    // The delete command uses singletons so we expect that directory to be used
    lazy val workflowLogPath = WorkflowLogger.workflowLogConfiguration.get.dir.resolve(s"workflow.$workflowId.log")
    copyWorkflowLogsActor ! CopyWorkflowLogsActor.Copy(workflowId, destinationPath)

    // Send a shutdown before the delete
    deathWatch.watch(copyWorkflowLogsActor)
    copyWorkflowLogsActor ! GracefulShutdownHelper.ShutdownCommand

    // Because the copy failed we instead expect a delete command to be sent to the ioActor
    val deleteCommand = DefaultIoDeleteCommand(workflowLogPath, swallowIOExceptions = true)
    ioActor.expectMsg(msgWait, deleteCommand)

    // Test that the actor is still alive and receiving messages even after a shutdown was requested
    EventFilter.error(pattern = "Failed to delete workflow logs", occurrences = 1).intercept {
      copyWorkflowLogsActor ! IoFailure(deleteCommand, new Exception("everything's fine, I am an expected delete fail") with NoStackTrace)
    }

    // Then the actor should shutdown
    deathWatch.expectTerminated(copyWorkflowLogsActor, msgWait)
  }

  it should "exit if creating a delete commands fails after failing to create a copy command" in {
    val workflowId = WorkflowId.randomId()
    val destinationPath = DefaultPathBuilder.createTempFile(s"test_file_$workflowId.", ".file", Option(tempDir))
    val partialIoCommandBuilder = new PartialIoCommandBuilder {
      override def copyCommand: PartialFunction[(Path, Path), Try[IoCopyCommand]] = {
        case _ => Failure(new Exception("everything's fine, I am an expected copy fail") with NoStackTrace)
      }

      override def deleteCommand: PartialFunction[(Path, Boolean), Try[IoDeleteCommand]] = {
        case _ => Failure(new Exception("everything's fine, I am an expected delete fail") with NoStackTrace)
      }
    }
    val ioCommandBuilder = new IoCommandBuilder(List(partialIoCommandBuilder))
    val props =
      CopyWorkflowLogsActor.props(
        serviceRegistryActor = serviceRegistryActor.ref,
        ioActor = ioActor.ref,
        workflowLogConfigurationOption = WorkflowLogger.workflowLogConfiguration,
        copyCommandBuilder = ioCommandBuilder,
        deleteCommandBuilder = ioCommandBuilder,
      )
    val copyWorkflowLogsActor = system.actorOf(props, "testCopyWorkflowLogsActorFailDelete")

    EventFilter.error(pattern = "Failed to delete workflow logs", occurrences = 1).intercept {
      EventFilter.error(pattern = "Failed to copy workflow logs", occurrences = 1).intercept {
        copyWorkflowLogsActor ! CopyWorkflowLogsActor.Copy(workflowId, destinationPath)
      }
    }

    // Send a shutdown
    deathWatch.watch(copyWorkflowLogsActor)
    copyWorkflowLogsActor ! GracefulShutdownHelper.ShutdownCommand

    // Then the actor should shutdown
    deathWatch.expectTerminated(copyWorkflowLogsActor, msgWait)
  }
}
