package cromwell.engine.finalcall

import cromwell.engine._
import cromwell.engine.backend.OldStyleWorkflowDescriptor
import cromwell.webservice.WorkflowMetadataResponse

import scala.concurrent.ExecutionContext
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object CopyWorkflowLogCall extends FinalCallCompanion[CopyWorkflowLogCall] {
  override val finalCallName = "copy_workflow_log"

  override def createCall(workflow: OldStyleWorkflowDescriptor) = CopyWorkflowLogCall(workflow)
}

/**
  * Final call implementation that copies the workflow log to a specified destination.
  */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class CopyWorkflowLogCall(override val workflow: OldStyleWorkflowDescriptor) extends OldStyleFinalCall {
  override val companion = CopyWorkflowLogCall
  override val handle = CopyWorkflowLogHandle

  override def execute(workflowMetadataResponse: WorkflowMetadataResponse)(implicit ec: ExecutionContext) = {
    workflow.copyWorkflowLog()
  }
}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case object CopyWorkflowLogHandle extends FinalCallHandle
