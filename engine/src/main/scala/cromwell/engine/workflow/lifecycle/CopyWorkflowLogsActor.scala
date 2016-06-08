package cromwell.engine.workflow.lifecycle

import java.nio.file.Path

import akka.actor.Props
import cromwell.core.{ExecutionStore, OutputStore, PathCopier, WorkflowId}
import cromwell.engine.EngineWorkflowDescriptor

object CopyWorkflowLogsActor {
  def props(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor, executionStore: ExecutionStore, outputStore: OutputStore) = Props(
    new CopyWorkflowLogsActor(workflowId, workflowDescriptor, executionStore, outputStore)
  )
}

class CopyWorkflowLogsActor(workflowId: WorkflowId, val workflowDescriptor: EngineWorkflowDescriptor,
                            executionStore: ExecutionStore, outputStore: OutputStore)
  extends EngineWorkflowCopyFinalizationActor {

  override def copyFiles(): Unit = {
    for {
      workflowLogDirString <- getWorkflowOption("workflow_log_dir")
      tempLogFilePath <- getTempLogFilePathOption
    } yield copyWorkflowLog(tempLogFilePath, workflowLogDirString)
  }

  // TODO: PBE: https://github.com/broadinstitute/cromwell/issues/814
  private def getTempLogFilePathOption: Option[Path] = ???

  private def copyWorkflowLog(tempLogFilePath: Path, workflowLogDirString: String): Unit = {
    val workflowLogDirPath = convertStringToPath(workflowLogDirString)
    val destinationFilePath = workflowLogDirPath.resolve(tempLogFilePath.getFileName.toString)
    PathCopier.copy(tempLogFilePath, destinationFilePath)
  }
}
