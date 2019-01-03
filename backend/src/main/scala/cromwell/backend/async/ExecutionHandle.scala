package cromwell.backend.async

import cromwell.backend.BackendJobDescriptor
import cromwell.backend.async.AsyncBackendJobExecutionActor.JobId
import cromwell.core.path.Path
import cromwell.core.{CallOutputs, ExecutionEvent}

/**
 * Trait to encapsulate whether an execution is complete and if so provide a result.  Useful in conjunction
 * with the `poll` API to feed results of previous job status queries forward.
 */
sealed trait ExecutionHandle {
  def isDone: Boolean
  def result: ExecutionResult
}

final case class PendingExecutionHandle[BackendJobId <: JobId, BackendRunInfo, BackendRunState]
(
  jobDescriptor: BackendJobDescriptor,
  pendingJob: BackendJobId,
  runInfo: Option[BackendRunInfo],
  previousState: Option[BackendRunState]
) extends ExecutionHandle {
  override val isDone = false
  override val result = NonRetryableExecution(new IllegalStateException("PendingExecutionHandle cannot yield a result"))
}

final case class SuccessfulExecutionHandle(outputs: CallOutputs, returnCode: Int, jobDetritusFiles: Map[String, Path], executionEvents: Seq[ExecutionEvent], resultsClonedFrom: Option[BackendJobDescriptor] = None) extends ExecutionHandle {
  override val isDone = true
  override val result = SuccessfulExecution(outputs, returnCode, jobDetritusFiles, executionEvents, resultsClonedFrom)
}

final case class FailedNonRetryableExecutionHandle(throwable: Throwable, returnCode: Option[Int] = None) extends ExecutionHandle {
  override val isDone = true
  override val result = NonRetryableExecution(throwable, returnCode)
}

final case class FailedRetryableExecutionHandle(throwable: Throwable, returnCode: Option[Int] = None) extends ExecutionHandle {
  override val isDone = true
  override val result = RetryableExecution(throwable, returnCode)
}

case object AbortedExecutionHandle extends ExecutionHandle {
  override def isDone: Boolean = true
  override def result: ExecutionResult = AbortedExecution
}
