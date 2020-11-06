package cromwell.engine.workflow.lifecycle.finalization

import java.io.IOException

import akka.actor.SupervisorStrategy.Restart
import akka.actor.{Actor, ActorLogging, ActorRef, OneForOneStrategy, Props}
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core._
import cromwell.core.io._
import cromwell.core.logging.WorkflowLogger
import cromwell.core.logging.WorkflowLogger.WorkflowLogConfiguration
import cromwell.core.path.Path
import cromwell.engine.workflow.WorkflowMetadataHelper
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success, Try}

object CopyWorkflowLogsActor {
  // Commands
  case class Copy(workflowId: WorkflowId, destinationDirPath: Path)

  val strategy: OneForOneStrategy = OneForOneStrategy(maxNrOfRetries = 3) {
    case _: IOException => Restart
  }

  def props(serviceRegistryActor: ActorRef,
            ioActor: ActorRef,
            workflowLogConfigurationOption: Option[WorkflowLogConfiguration] =
            WorkflowLogger.workflowLogConfiguration,
            /*
            The theory is that the `GcsBatchCommandBuilder` copies the temporary workflow logs from the local disk to
            GCS. Then later, the separate `DefaultIOCommandBuilder` deletes files from the local disk.

            But...

            I believe the `GcsBatchCommandBuilder` _always_ fails to create a copy command as it only copies from GCS to
            GCS. Since a GCS to GCS copy command isn't created the `IoCommandBuilder.copyCommand` returns its
            `DefaultIoCopyCommand`.

            Then `cromwell.engine.io.nio.NioFlow.copy` does the copy.

            So in theory one could just use the `DefaultIoCommandBuilder` for both command builders. However
            attempting-to-use the `GcsBatchCommandBuilder` to create copy commands is the way this was previously
            implemented, It Works (TM), and I'm not changing it for now.
             */
            copyCommandBuilder: IoCommandBuilder = GcsBatchCommandBuilder,
            deleteCommandBuilder: IoCommandBuilder = DefaultIoCommandBuilder,
           ): Props = {
    Props(new CopyWorkflowLogsActor(
      serviceRegistryActor = serviceRegistryActor,
      ioActor = ioActor,
      workflowLogConfigurationOption = workflowLogConfigurationOption,
      copyCommandBuilder = copyCommandBuilder,
      deleteCommandBuilder = deleteCommandBuilder,
    )).withDispatcher(IoDispatcher)
  }
}

// This could potentially be turned into a more generic "Copy/Move something from A to B"
// Which could be used for other copying work (outputs, call logs..)
class CopyWorkflowLogsActor(override val serviceRegistryActor: ActorRef,
                            override val ioActor: ActorRef,
                            workflowLogConfigurationOption: Option[WorkflowLogConfiguration],
                            copyCommandBuilder: IoCommandBuilder,
                            deleteCommandBuilder: IoCommandBuilder,
                           ) extends Actor
  with ActorLogging with IoClientHelper with WorkflowMetadataHelper with MonitoringCompanionHelper {

  implicit val ec: ExecutionContext = context.dispatcher

  private def tryCopyLog(src: Path, dest: Path, workflowId: WorkflowId): Try[Unit] = {
    // In order to keep "copy and then delete" operations atomic as far as monitoring is concerned, removeWork will only be called
    // when the delete is complete (successfully or not), or when the copy completes if WorkflowLogger.isTemporary is false
    addWork()
    // Send the workflowId as context along with the copy so we can update metadata when the response comes back
    copyCommandBuilder.copyCommand(src, dest).map(sendIoCommandWithContext(_, workflowId))
  }

  private def deleteLog(src: Path): Unit = if (workflowLogConfigurationOption.exists(_.temporary)) {
    deleteCommandBuilder.deleteCommand(src) match {
      case Success(command) => sendIoCommand(command)
      case Failure(failure) =>
        log.error(s"Failed to delete workflow logs from ${src.pathAsString}: ${failure.getMessage}")
        removeWork()
    }
  } else removeWork()
  
  private def updateLogsPathInMetadata(workflowId: WorkflowId, path: Path): Unit = {
    val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.WorkflowLog), MetadataValue(path.pathAsString))
    serviceRegistryActor ! PutMetadataAction(metadataEventMsg)
  }

  private def copyLogsReceive: Receive = {
    case CopyWorkflowLogsActor.Copy(workflowId, destinationDir) =>
      val workflowLogger = new WorkflowLogger(
        loggerName = self.path.name,
        workflowId = workflowId.toPossiblyNotRoot,
        rootWorkflowId = workflowId.toRoot,
        akkaLogger = Option(log)
      )

      workflowLogger.workflowLogPath foreach { src =>
        if (Try(src.exists).getOrElse(false)) {
          val destPath = destinationDir.resolve(src.name)
          workflowLogger.info(s"Copying workflow logs from $src to $destPath")

          tryCopyLog(src, destPath, workflowId) match {
            case Failure(failure) =>
              log.error(
                cause = failure,
                message =
                  s"Failed to copy workflow logs from ${src.pathAsString} to ${destPath.pathAsString}: " +
                    s"${failure.getMessage}",
              )
              deleteLog(src)
            case Success(_) =>
              // Deliberately not deleting the file here, that will be done in batch in `deleteLog`
              // after the copy is terminal.
          }
          workflowLogger.close()
        }
      }
      
    case (workflowId: WorkflowId, IoSuccess(copy: IoCopyCommand, _)) =>
      updateLogsPathInMetadata(workflowId, copy.destination)
      deleteLog(copy.source)
      
    case (workflowId: WorkflowId, IoFailAck(copy: IoCopyCommand, failure)) =>
      pushWorkflowFailures(workflowId, List(new IOException("Could not copy workflow logs", failure)))
      log.error(failure, s"Failed to copy workflow logs from ${copy.source.pathAsString} to ${copy.destination.pathAsString}")
      deleteLog(copy.source)
      
    case IoSuccess(_: IoDeleteCommand, _) => removeWork()
      
    case IoFailAck(delete: IoDeleteCommand, failure) =>
      removeWork()
      log.error(failure, s"Failed to delete workflow logs from ${delete.file.pathAsString}")

    case other => log.warning(s"CopyWorkflowLogsActor received an unexpected message: $other")
  }

  /*_*/
  // shh intellij is ok...
  // https://stackoverflow.com/questions/36679973/controlling-false-intellij-code-editor-error-in-scala-plugin
  override def receive: Receive = monitoringReceive orElse ioReceive orElse copyLogsReceive
  /*_*/

  override def preRestart(t: Throwable, message: Option[Any]): Unit = {
    message foreach self.forward
  }

  override protected def onTimeout(message: Any, to: ActorRef): Unit = message match {
    case copy: IoCopyCommand =>
      log.error(s"Failed to copy workflow logs from ${copy.source.pathAsString} to ${copy.destination.pathAsString}: Timeout")
      deleteLog(copy.source)
    case delete: IoDeleteCommand =>
      log.error(s"Failed to delete workflow logs from ${delete.file.pathAsString}: Timeout")
    case _ =>
  }
}
