package cromwell.engine.finalcall

import cromwell.engine.backend.WorkflowDescriptor
import cromwell.webservice.WorkflowMetadataResponse

/**
  * Final call implementation that copies the workflow log to a specified destination.
  */
object CopyWorkflowLogCall extends CopyingFinalCall {
  override val finalCallName = "copy_workflow_log"

  override def createFinalCallCopies(workflow: WorkflowDescriptor, metadata: WorkflowMetadataResponse) = {
    val workflowLogPath = workflow.workflowLogPath
    val workflowLogDir = FinalCall.getWorkflowOption(workflow, WorkflowDescriptor.WorkflowLogDirOptionKey)
    val copyOption = for {
      src <- workflowLogPath
      dir <- workflowLogDir
      dest = workflow.toFileSystemPath(dir).resolve(workflow.workflowLogName)
    } yield FinalCallCopy(src, dest)
    copyOption.toSeq
  }
}
