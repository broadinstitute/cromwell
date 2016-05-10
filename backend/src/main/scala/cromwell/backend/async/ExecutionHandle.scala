package cromwell.backend.async

import cromwell.backend.{BackendJobDescriptor, ExecutionHash}
import cromwell.core.CallOutputs

/**
 * Trait to encapsulate whether an execution is complete and if so provide a result.  Useful in conjunction
 * with the `poll` API to feed results of previous job status queries forward.
 */
trait ExecutionHandle {
  def isDone: Boolean
  def result: ExecutionResult
}

final case class SuccessfulExecutionHandle(outputs: CallOutputs, returnCode: Int, hash: ExecutionHash, resultsClonedFrom: Option[BackendJobDescriptor] = None) extends ExecutionHandle {
  override val isDone = true
  override val result = SuccessfulExecution(outputs, returnCode, hash, resultsClonedFrom)
}

final case class FailedNonRetryableExecutionHandle(throwable: Throwable, returnCode: Option[Int] = None) extends ExecutionHandle {
  override val isDone = true
  override val result = new NonRetryableExecution(throwable, returnCode)
}

final case class FailedRetryableExecutionHandle(throwable: Throwable, returnCode: Option[Int] = None) extends ExecutionHandle {
  override val isDone = true
  override val result = new RetryableExecution(throwable, returnCode)
}

case object AbortedExecutionHandle extends ExecutionHandle {
  override def isDone: Boolean = true
  override def result: ExecutionResult = AbortedExecution
}
