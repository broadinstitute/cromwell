package cromwell.backend.async

import cromwell.backend.BackendJobDescriptor
import cromwell.core.JobOutputs

/**
 * ADT representing the result of an execution of a BackendCall.
 */
sealed trait ExecutionResult

/**
 * A successful execution with resolved outputs.
 */
final case class SuccessfulExecution(outputs: JobOutputs, returnCode: Int, jobOutputFiles: Map[String, String], resultsClonedFrom: Option[BackendJobDescriptor] = None) extends ExecutionResult

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
final case class NonRetryableExecution(e: Throwable, returnCode: Option[Int] = None) extends FailedExecution
final case class RetryableExecution(e: Throwable, returnCode: Option[Int] = None) extends FailedExecution
