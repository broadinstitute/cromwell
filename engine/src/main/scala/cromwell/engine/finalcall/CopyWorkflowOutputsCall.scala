package cromwell.engine.finalcall

import cromwell.engine.backend.OldStyleWorkflowDescriptor
import cromwell.webservice.WorkflowMetadataResponse

import scala.concurrent.ExecutionContext
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object CopyWorkflowOutputsCall extends FinalCallCompanion[CopyWorkflowOutputsCall] {
  override val finalCallName = "copy_workflow_outputs"

  override def createCall(workflow: OldStyleWorkflowDescriptor) = CopyWorkflowOutputsCall(workflow)
}

/**
  * Final call implementation that copies workflow outputs to a specified destination.
  */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class CopyWorkflowOutputsCall(override val workflow: OldStyleWorkflowDescriptor) extends OldStyleFinalCall {
  override val companion = CopyWorkflowOutputsCall
  override val handle = CopyWorkflowOutputsHandle

  override def execute(workflowMetadataResponse: WorkflowMetadataResponse)(implicit ec: ExecutionContext) = {
    workflow.copyWorkflowOutputs(workflowMetadataResponse)
  }
}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case object CopyWorkflowOutputsHandle extends FinalCallHandle
