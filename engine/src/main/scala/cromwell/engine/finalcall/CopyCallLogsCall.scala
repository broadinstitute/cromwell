package cromwell.engine.finalcall

import cromwell.engine.backend.OldStyleWorkflowDescriptor
import cromwell.webservice.WorkflowMetadataResponse

import scala.concurrent.ExecutionContext
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object CopyCallLogsCall extends FinalCallCompanion[CopyCallLogsCall] {
  override val finalCallName = "copy_call_logs"

  override def createCall(workflow: OldStyleWorkflowDescriptor) = CopyCallLogsCall(workflow)
}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class CopyCallLogsCall(override val workflow: OldStyleWorkflowDescriptor) extends OldStyleFinalCall {
  override val companion = CopyCallLogsCall
  override val handle = CopyCallLogsCallHandle

  override def execute(workflowMetadataResponse: WorkflowMetadataResponse)(implicit ec: ExecutionContext) = {
    workflow.copyCallLogs(workflowMetadataResponse)
  }
}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case object CopyCallLogsCallHandle extends FinalCallHandle
