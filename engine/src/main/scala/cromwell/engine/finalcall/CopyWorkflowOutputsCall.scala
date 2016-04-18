package cromwell.engine.finalcall

import cromwell.engine.backend.WorkflowDescriptor
import cromwell.webservice.WorkflowMetadataResponse

/**
  * Final call implementation that copies workflow outputs to a specified destination.
  */
object CopyWorkflowOutputsCall extends CopyingFinalCall {
  override val finalCallName = "copy_workflow_outputs"

  override def createFinalCallCopies(workflow: WorkflowDescriptor, metadata: WorkflowMetadataResponse) = {
    val workflowOutputsDir = FinalCall.getWorkflowOption(workflow, WorkflowDescriptor.WorkflowOutputsOptionKey)
    workflowOutputsDir.toSeq flatMap { dir =>
      val outputs = metadata.outputs.toIndexedSeq.flatMap(_.values)
      outputs flatMap CopyingFinalCall.copyToDir(workflow, dir)
    }
  }
}
