package cromwell.engine.backend

import cromwell.backend.ExecutionHash
import cromwell.core.JobOutputs

/**
 * Trait to encapsulate whether an execution is complete and if so provide a result.  Useful in conjunction
 * with the `poll` API to feed results of previous job status queries forward.
 */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
trait OldStyleExecutionHandle {
  def isDone: Boolean
  def result: OldStyleExecutionResult
}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
final case class CompletedExecutionHandle(override val result: OldStyleExecutionResult) extends OldStyleExecutionHandle {
  override val isDone = true
}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
final case class SuccessfulExecutionHandle(outputs: JobOutputs, returnCode: Int, hash: ExecutionHash, resultsClonedFrom: Option[OldStyleBackendCallJobDescriptor] = None) extends OldStyleExecutionHandle {
  override val isDone = true
  override val result = OldStyleSuccessfulBackendCallExecution(outputs, returnCode, hash, resultsClonedFrom)
}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
final case class FailedExecutionHandle(throwable: Throwable, returnCode: Option[Int] = None) extends OldStyleExecutionHandle {
  override val isDone = true
  override val result = new OldStyleNonRetryableFailedExecution(throwable, returnCode)
}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
final case class RetryableExecutionHandle(throwable: Throwable, returnCode: Option[Int] = None) extends OldStyleExecutionHandle {
  override val isDone = true
  override val result = new OldStyleRetryableFailedExecution(throwable, returnCode)
}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case object AbortedExecutionHandle extends OldStyleExecutionHandle {
  override def isDone: Boolean = true
  override def result: OldStyleExecutionResult = OldStyleAbortedExecution
}
