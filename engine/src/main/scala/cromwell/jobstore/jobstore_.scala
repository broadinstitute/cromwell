package cromwell.jobstore

import cromwell.backend.BackendJobDescriptorKey
import cromwell.backend.BackendJobExecutionActor.{
  FetchedFromJobStore,
  JobFailedNonRetryableResponse,
  JobFailedRetryableResponse,
  JobSucceededResponse
}
import cromwell.core.{CallOutputs, WorkflowId}

case class JobStoreKey(workflowId: WorkflowId, callFqn: String, index: Option[Int], attempt: Int) {
  private lazy val indexString = index map { _.toString } getOrElse "NA"
  lazy val tag = s"$workflowId:$callFqn:$indexString:$attempt"
}

sealed trait JobResult {
  def toBackendJobResponse(key: BackendJobDescriptorKey) = this match {
    // Always puts `None` for `dockerImageUsed` for a successfully completed job on restart. This shouldn't be a
    // problem since `saveJobCompletionToJobStore` in EJEA will already have sent this to metadata.
    case JobResultSuccess(returnCode, jobOutputs) =>
      JobSucceededResponse(key,
                           returnCode,
                           jobOutputs,
                           None,
                           Seq.empty,
                           None,
                           resultGenerationMode = FetchedFromJobStore
      )
    case JobResultFailure(returnCode, reason, false) => JobFailedNonRetryableResponse(key, reason, returnCode)
    case JobResultFailure(returnCode, reason, true) => JobFailedRetryableResponse(key, reason, returnCode)
  }
}
case class JobResultSuccess(returnCode: Option[Int], jobOutputs: CallOutputs) extends JobResult
case class JobResultFailure(returnCode: Option[Int], reason: Throwable, retryable: Boolean) extends JobResult
