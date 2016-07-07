package cromwell

import cromwell.core.{JobKey, JobOutputs, WorkflowId}

/**
  * Created by chrisl on 7/5/16.
  */
package object jobstore {
  case class JobStoreKey(workflowId: WorkflowId, callFqn: String, index: Option[Int], attempt: Int)

  sealed trait JobResult
  case class JobResultSuccess(returnCode: Option[Int], jobOutputs: JobOutputs) extends JobResult
  case class JobResultFailure(returnCode: Option[Int], reason: Throwable) extends JobResult

}
