package cromwell.engine.finalcall

import cromwell.engine._
import cromwell.webservice.WorkflowMetadataResponse

import scala.concurrent.ExecutionContext

object CopyWorkflowOutputsCall extends FinalCallCompanion[CopyWorkflowOutputsCall] {
  override val finalCallName = "copy_workflow_outputs"

  override def createCall(workflow: WorkflowDescriptor) = CopyWorkflowOutputsCall(workflow)
}

/**
  * Final call implementation that copies workflow outputs to a specified destination.
  */
case class CopyWorkflowOutputsCall(override val workflow: WorkflowDescriptor) extends FinalCall {
  override val companion = CopyWorkflowOutputsCall
  override val handle = CopyWorkflowOutputsHandle

  override def execute(workflowMetadataResponse: WorkflowMetadataResponse)(implicit ec: ExecutionContext) = {
    workflow.copyWorkflowOutputs(workflowMetadataResponse)
  }
}

case object CopyWorkflowOutputsHandle extends FinalCallHandle
