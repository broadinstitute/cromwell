package cromwell.engine.backend

import cromwell.engine.ExecutionEventEntry
import cromwell.engine.ExecutionHash
import cromwell.engine.CallOutputs

/**
 * ADT representing the result of an execution of a BackendCall.
 */
sealed trait ExecutionResult

case object SuccessfulFinalCallExecution extends ExecutionResult

/**
 * A successful execution with resolved outputs.
 */
final case class SuccessfulBackendCallExecution(outputs: CallOutputs, executionEvents: Seq[ExecutionEventEntry], returnCode: ScriptReturnCode, hash: ExecutionHash, resultsClonedFrom: Option[BackendCall] = None) extends ExecutionResult

/**
 * A user-requested abort of the command.
 */
case object AbortedExecution extends ExecutionResult

sealed trait FailedExecution extends ExecutionResult {
  def e: Throwable
  def returnCode: Option[ReturnCode]
}

/**
  * Failed execution, possibly having a return code.
  */
final case class NonRetryableExecution(e: Throwable, returnCode: Option[ReturnCode] = None, events: Seq[ExecutionEventEntry] = Seq.empty) extends FailedExecution
final case class RetryableExecution(e: Throwable, returnCode: Option[ReturnCode] = None, events: Seq[ExecutionEventEntry] = Seq.empty) extends FailedExecution
