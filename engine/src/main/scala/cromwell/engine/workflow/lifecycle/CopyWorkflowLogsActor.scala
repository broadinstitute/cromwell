package cromwell.engine.workflow.lifecycle

import java.nio.file.Path

import akka.actor.{Actor, ActorLogging, Props}
import better.files._
import cromwell.core._
import cromwell.core.logging.WorkflowLogger
import cromwell.database.obj.WorkflowMetadataKeys
import cromwell.services.MetadataServiceActor.PutMetadataAction
import cromwell.services.{MetadataEvent, MetadataKey, MetadataValue, ServiceRegistryClient}

object CopyWorkflowLogsActor {
  // Commands
  case class Copy(workflowId: WorkflowId, destinationDirPath: Path)

  def props = Props(new CopyWorkflowLogsActor())
}

// This could potentially be turned into a more generic "Copy/Move something from A to B"
// Which could be used for other copying work (outputs, call logs..)
class CopyWorkflowLogsActor extends Actor with ActorLogging with PathFactory with ServiceRegistryClient {

  // Return the destination Path if copying happened
  def copyAndClean(dest: Path, workflowLogger: WorkflowLogger): Option[Path] = {
    workflowLogger.workflowLogPath map { src =>
      dest.createDirectories()

      val destPath = dest.resolve(src.getFileName)
      workflowLogger.info(s"Copying workflow logs from ${src.toAbsolutePath} to $destPath")

      src.copyTo(destPath)
      if (WorkflowLogger.isTemporary) src.delete()

      destPath
    }
  }

  override def receive = {
    case CopyWorkflowLogsActor.Copy(workflowId, destinationDir) =>

      val workflowLogger = new WorkflowLogger(self.path.name, workflowId, Option(log))

      copyAndClean(destinationDir, workflowLogger) foreach { destinationPath =>
        val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.WorkflowLog), MetadataValue(destinationPath))
        serviceRegistryActor ! PutMetadataAction(metadataEventMsg)
      }

      context stop self
  }
}
