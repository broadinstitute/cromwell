package cromwell.core.abort

import cromwell.core.WorkflowId

sealed trait AbortResponse {
  val workflowId: WorkflowId
}

final case class WorkflowAbortFailureResponse(workflowId: WorkflowId, failure: Throwable) extends AbortResponse

sealed trait SuccessfulAbortResponse extends AbortResponse
final case class WorkflowAbortedResponse(workflowId: WorkflowId) extends SuccessfulAbortResponse
final case class WorkflowAbortRequestedResponse(workflowId: WorkflowId) extends SuccessfulAbortResponse
