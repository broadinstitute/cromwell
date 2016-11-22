package cromwell.jobstore

import cromwell.core.{WorkflowId, _}

case class JobStoreKey(workflowId: WorkflowId, callFqn: String, index: Option[Int], attempt: Int)

sealed trait JobResult
case class JobResultSuccess(returnCode: Option[Int], jobOutputs: CallOutputs) extends JobResult
case class JobResultFailure(returnCode: Option[Int], reason: Throwable, retryable: Boolean) extends JobResult

