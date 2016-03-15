package cromwell.engine.backend

import cromwell.engine.{CallOutputs, ExecutionEventEntry, ExecutionHash}

/**
 * Trait to encapsulate whether an execution is complete and if so provide a result.  Useful in conjunction
 * with the `poll` API to feed results of previous job status queries forward.
 */
trait ExecutionHandle {
  def isDone: Boolean
  def result: ExecutionResult
}

final case class CompletedExecutionHandle(override val result: ExecutionResult) extends ExecutionHandle {
  override val isDone = true
}

final case class SuccessfulExecutionHandle(outputs: CallOutputs, events: Seq[ExecutionEventEntry], returnCode: Int, hash: ExecutionHash, resultsClonedFrom: Option[BackendCallJobDescriptor] = None) extends ExecutionHandle {
  override val isDone = true
  override val result = SuccessfulBackendCallExecution(outputs, events, returnCode, hash, resultsClonedFrom)
}

final case class FailedExecutionHandle(throwable: Throwable, returnCode: Option[Int] = None, events: Seq[ExecutionEventEntry] = Seq.empty) extends ExecutionHandle {
  override val isDone = true
  override val result = new NonRetryableExecution(throwable, returnCode, events)
}

final case class RetryableExecutionHandle(throwable: Throwable, returnCode: Option[Int] = None, events: Seq[ExecutionEventEntry] = Seq.empty) extends ExecutionHandle {
  override val isDone = true
  override val result = new RetryableExecution(throwable, returnCode, events)
}

case object AbortedExecutionHandle extends ExecutionHandle {
  override def isDone: Boolean = true
  override def result: ExecutionResult = AbortedExecution
}