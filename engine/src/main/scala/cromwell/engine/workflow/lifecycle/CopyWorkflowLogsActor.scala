package cromwell.engine.workflow.lifecycle

import java.nio.file.Path

import akka.actor.{Actor, Props}
import better.files._
import cromwell.core.WorkflowOptions._
import cromwell.core._
import cromwell.core.logging.{WorkflowLogger, WorkflowLogging}
import cromwell.database.obj.WorkflowMetadataKeys
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.services.MetadataServiceActor.PutMetadataAction
import cromwell.services.{MetadataValue, MetadataKey, MetadataEvent, ServiceRegistryClient}

object CopyWorkflowLogsActor {
  // Commands
  case object Start

  def props(workflowDescriptor: EngineWorkflowDescriptor) = Props(
    new CopyWorkflowLogsActor(workflowDescriptor)
  )
}

class CopyWorkflowLogsActor(val workflowDescriptor: EngineWorkflowDescriptor) extends Actor with PathFactory with WorkflowLogging with ServiceRegistryClient {

  override lazy val workflowId = workflowDescriptor.id
  val destinationPath = {
    workflowDescriptor.getWorkflowOption(FinalWorkflowLogDir) map { path =>
      buildPath(path, workflowDescriptor.engineFilesystems)
    }
  }

  def copyAndClean(dest: Path) = {
    workflowLogger.workflowLogPath foreach { src =>
      val destPath = dest.toAbsolutePath.resolve(src.getFileName)
      workflowLogger.info(s"Copying workflow logs from ${src.toAbsolutePath} to $destPath")

      PathCopier.copy(src.toAbsolutePath, destPath)
      if (WorkflowLogger.workflowLogConfiguration exists { _.temporary }) src.delete(ignoreIOExceptions = false)

      val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.WorkflowLog), MetadataValue(destPath.toAbsolutePath))
      serviceRegistryActor ! PutMetadataAction(metadataEventMsg)
    }
  }

  override def receive = {
    case CopyWorkflowLogsActor.Start =>
      destinationPath foreach copyAndClean
      context stop self
  }
}
