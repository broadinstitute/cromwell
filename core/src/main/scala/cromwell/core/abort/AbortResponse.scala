package cromwell.core.abort

import cromwell.core.WorkflowId

sealed trait AbortResponse
case class WorkflowAbortingResponse(workflowId: WorkflowId, restarted: Boolean) extends AbortResponse
case class WorkflowAbortFailureResponse(workflowId: WorkflowId, failure: Throwable) extends AbortResponse
