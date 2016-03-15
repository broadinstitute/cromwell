package cromwell.engine.backend

import cromwell.engine.{CallOutputs, ExecutionEventEntry, ExecutionHash}

/**
 * Represents a Call that is intended to be run on a specific Backend.
 *
 * A BackendCall is created by calling `Backend.bindCall(c: Call)`.  The
 * combination of a Call and a Backend give enough context to be able to
 * instantiate the Call's Task's command line, and provide an implementation
 * of the WDL standard library functions as well as a lookup functions to
 * resolve identifiers in expressions.
 *
 * The BackendCall also includes a Map[String, WdlValue] (CallInputs) which
 * represents the locally-qualified input values for this call.
 *
 * BackendCall aims to package everything that a Call needs to run into one
 * object so a client can do a parameterless `.execute()` to run the job.
 *
 * To understand why one needs a Backend and locally-qualified inputs in order
 * to instantiate a command, consider this example:
 *
 * task test {
 *   Array[String] something
 *   command {
 *     ./analysis --something-list=${write_lines(something)}
 *   }
 * }
 *
 * Without the Backend and the locallyQualifiedInputs, the expression
 * `write_lines(something)` could not be evaluated because:
 *
 * 1) We need to resolve "something" in locallyQualifiedInputs
 * 2) We need a place to write the results of `write_lines()`, and that location
 *    is backend-specific
 */

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