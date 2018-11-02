package cromwell.core.abort

import cromwell.core.WorkflowId

sealed trait AbortResponse {
  val workflowId: WorkflowId
}

final case class WorkflowAbortFailureResponse(workflowId: WorkflowId, failure: Throwable) extends AbortResponse

sealed trait SuccesfulAbortResponse extends AbortResponse
final case class WorkflowAbortedResponse(workflowId: WorkflowId) extends SuccesfulAbortResponse
final case class WorkflowAbortRequestedResponse(workflowId: WorkflowId) extends SuccesfulAbortResponse
