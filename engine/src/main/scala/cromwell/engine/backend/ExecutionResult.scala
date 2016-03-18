package cromwell.engine.backend

import cromwell.core.CallOutputs
import cromwell.engine.ExecutionEventEntry

/**
 * ADT representing the result of an execution of a BackendCall.
 */
sealed trait ExecutionResult

case object SuccessfulFinalCallExecution extends ExecutionResult

/**
 * A successful execution with resolved outputs.
 */
final case class SuccessfulBackendCallExecution(outputs: CallOutputs, executionEvents: Seq[ExecutionEventEntry], returnCode: Int, hash: ExecutionHash, resultsClonedFrom: Option[BackendCallJobDescriptor] = None) extends ExecutionResult

/**
 * A user-requested abort of the command.
 */
case object AbortedExecution extends ExecutionResult

sealed trait FailedExecution extends ExecutionResult {
  def e: Throwable
  def returnCode: Option[Int]
}

/**
  * Failed execution, possibly having a return code.
  */
final case class NonRetryableExecution(e: Throwable, returnCode: Option[Int] = None, events: Seq[ExecutionEventEntry] = Seq.empty) extends FailedExecution
final case class RetryableExecution(e: Throwable, returnCode: Option[Int] = None, events: Seq[ExecutionEventEntry] = Seq.empty) extends FailedExecution
