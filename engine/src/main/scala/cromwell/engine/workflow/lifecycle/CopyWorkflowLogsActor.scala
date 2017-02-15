package cromwell.engine.workflow.lifecycle

import java.io.IOException

import akka.actor.SupervisorStrategy.Restart
import akka.actor.{Actor, ActorLogging, ActorRef, OneForOneStrategy, Props}
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core._
import cromwell.core.actor.RobustClientHelper
import cromwell.core.io.{IoCopyCommand, IoDeleteCommand, IoFailure, IoSuccess}
import cromwell.core.logging.WorkflowLogger
import cromwell.core.path.Path
import cromwell.filesystems.gcs.GcsBatchCommandBuilder
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}

object CopyWorkflowLogsActor {
  // Commands
  case class Copy(workflowId: WorkflowId, destinationDirPath: Path)

  val strategy = OneForOneStrategy(maxNrOfRetries = 3) {
    case _: IOException => Restart
  }

  def props(serviceRegistryActor: ActorRef, ioActor: ActorRef) = Props(new CopyWorkflowLogsActor(serviceRegistryActor, ioActor)).withDispatcher(IoDispatcher)
}

// This could potentially be turned into a more generic "Copy/Move something from A to B"
// Which could be used for other copying work (outputs, call logs..)
class CopyWorkflowLogsActor(serviceRegistryActor: ActorRef, ioActor: ActorRef) extends Actor with ActorLogging with GcsBatchCommandBuilder with RobustClientHelper {

  def copyLog(src: Path, dest: Path, workflowId: WorkflowId) = {
    dest.parent.createPermissionedDirectories()
    // Send the workflowId as context along with the copy so we can update metadata when the response comes back
    ioActor ! (workflowId -> copyCommand(src, dest, overwrite = true))
  }

  def deleteLog(src: Path) = if (WorkflowLogger.isTemporary) {
    ioActor ! deleteCommand(src)
  }
  
  def updateLogsPathInMetadata(workflowId: WorkflowId, path: Path) = {
    val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.WorkflowLog), MetadataValue(path.pathAsString))
    serviceRegistryActor ! PutMetadataAction(metadataEventMsg)
  }

  def copyLogsReceive: Receive = {
    case CopyWorkflowLogsActor.Copy(workflowId, destinationDir) =>
      val workflowLogger = new WorkflowLogger(self.path.name, workflowId, Option(log))

      workflowLogger.workflowLogPath foreach { src =>
        if (src.exists) {
          val destPath = destinationDir.resolve(src.name)
          workflowLogger.info(s"Copying workflow logs from $src to $destPath")

          copyLog(src, destPath, workflowId)
        }
      }
      
    case (workflowId: WorkflowId, IoSuccess(copy: IoCopyCommand, _)) =>
      updateLogsPathInMetadata(workflowId, copy.destination)
      deleteLog(copy.source)
      
    case (workflowId, IoFailure(copy: IoCopyCommand, failure)) => 
      // TODO: do we still want to send something to metadata ? If so, what ?
      log.error(failure, s"Failed to copy workflow logs from ${copy.source.pathAsString} to ${copy.destination.pathAsString}")
      deleteLog(copy.source)
      
    case IoSuccess(_: IoDeleteCommand, _) => // Good !
      
    case IoFailure(delete: IoDeleteCommand, failure) =>
      log.error(failure, s"Failed to delete workflow logs from ${delete.file.pathAsString}")

    case other => log.warning(s"CopyWorkflowLogsActor received an unexpected message: $other")
  }
  
  override def receive = robustReceive orElse copyLogsReceive

  override def preRestart(t: Throwable, message: Option[Any]) = {
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
