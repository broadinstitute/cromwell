package cromwell

import cromwell.core.{JobKey, JobOutputs, WorkflowId}

package object jobstore {
  case class JobStoreKey(workflowId: WorkflowId, callFqn: String, index: Option[Int], attempt: Int)

  sealed trait JobResult
  case class JobResultSuccess(returnCode: Option[Int], jobOutputs: JobOutputs) extends JobResult
  case class JobResultFailure(returnCode: Option[Int], reason: Throwable, retryable: Boolean) extends JobResult

  implicit class EnhancedJobKey(val jobKey: JobKey) extends AnyVal {
    def toJobStoreKey(workflowId: WorkflowId): JobStoreKey = JobStoreKey(workflowId, jobKey.scope.fullyQualifiedName, jobKey.index, jobKey.attempt)
  }
}
