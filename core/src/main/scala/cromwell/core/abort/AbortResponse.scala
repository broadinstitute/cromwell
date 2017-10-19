package cromwell.core.abort

import cromwell.core.WorkflowId

sealed trait AbortResponse
case class WorkflowAbortingResponse(workflowId: WorkflowId) extends AbortResponse
case class WorkflowAbortFailureResponse(workflowId: WorkflowId, failure: Throwable) extends AbortResponse
