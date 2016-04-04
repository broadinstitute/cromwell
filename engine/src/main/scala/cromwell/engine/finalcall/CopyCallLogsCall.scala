package cromwell.engine.finalcall

import cromwell.engine.backend.WorkflowDescriptor
import cromwell.webservice.WorkflowMetadataResponse

object CopyCallLogsCall extends CopyingFinalCall {
  override val finalCallName = "copy_call_logs"

  override def createFinalCallCopies(workflow: WorkflowDescriptor, metadata: WorkflowMetadataResponse) = {
    val callLogsDir = FinalCall.getWorkflowOption(workflow, WorkflowDescriptor.CallLogsDirOptionKey)
    callLogsDir.toSeq flatMap copyCallLogs(workflow, metadata)
  }

  private def copyCallLogs(workflow: WorkflowDescriptor,
                            metadata: WorkflowMetadataResponse)(dir: String): Iterable[FinalCallCopy] = {
    val callLogs = for {
      callMetadatas <- metadata.calls.values
      callMetadata <- callMetadatas
      stdout = callMetadata.stdout.toSeq
      stderr = callMetadata.stderr.toSeq
      backend = callMetadata.backendLogs.toSeq.flatMap(_.values)
      callLog <- stdout ++ stderr ++ backend
    } yield callLog
    callLogs flatMap CopyingFinalCall.copyToDir(workflow, dir)
  }
}
