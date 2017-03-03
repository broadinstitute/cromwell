package cromwell.engine.workflow.lifecycle

import java.io.IOException

import akka.actor.SupervisorStrategy.Restart
import akka.actor.{Actor, ActorLogging, ActorRef, OneForOneStrategy, Props}
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core._
import cromwell.core.logging.WorkflowLogger
import cromwell.core.path.Path
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}

object CopyWorkflowLogsActor {
  // Commands
  case class Copy(workflowId: WorkflowId, destinationDirPath: Path)

  val strategy = OneForOneStrategy(maxNrOfRetries = 3) {
    case _: IOException => Restart
  }

  def props(serviceRegistryActor: ActorRef) = Props(new CopyWorkflowLogsActor(serviceRegistryActor)).withDispatcher(IoDispatcher)
}

// This could potentially be turned into a more generic "Copy/Move something from A to B"
// Which could be used for other copying work (outputs, call logs..)
class CopyWorkflowLogsActor(serviceRegistryActor: ActorRef)
    extends Actor
    with ActorLogging {

  def copyAndClean(src: Path, dest: Path) = {
    dest.parent.createPermissionedDirectories()

    src.copyTo(dest, overwrite = true)
    if (WorkflowLogger.isTemporary) {
      src.delete()
    }
  }

  override def receive = {
    case CopyWorkflowLogsActor.Copy(workflowId, destinationDir) =>
      val workflowLogger = new WorkflowLogger(self.path.name, workflowId, Option(log))

      workflowLogger.workflowLogPath foreach { src =>
        if (src.exists) {
          val destPath = destinationDir.resolve(src.name)
          workflowLogger.info(s"Copying workflow logs from $src to $destPath")

          copyAndClean(src, destPath)

          val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.WorkflowLog), MetadataValue(destPath))
          serviceRegistryActor ! PutMetadataAction(metadataEventMsg)
        }
      }
  }

  override def preRestart(t: Throwable, message: Option[Any]) = {
    message foreach self.forward
  }
}
