package cromwell.engine.finalcall

import cromwell.engine.WorkflowDescriptor
import cromwell.webservice.WorkflowMetadataResponse

import scala.concurrent.ExecutionContext

object CopyCallLogsCall extends FinalCallCompanion[CopyCallLogsCall] {
  override val finalCallName = "copy_call_logs"

  override def createCall(workflow: WorkflowDescriptor) = CopyCallLogsCall(workflow)
}

case class CopyCallLogsCall(override val workflow: WorkflowDescriptor) extends FinalCall {
  override val companion = CopyCallLogsCall
  override val handle = CopyCallLogsCallHandle

  override def execute(workflowMetadataResponse: WorkflowMetadataResponse)(implicit ec: ExecutionContext) = {
    workflow.copyCallLogs(workflowMetadataResponse)
  }
}

case object CopyCallLogsCallHandle extends FinalCallHandle
