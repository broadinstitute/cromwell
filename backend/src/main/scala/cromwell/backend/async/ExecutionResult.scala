package cromwell.backend.async

import java.nio.file.Path

import cromwell.backend.BackendJobDescriptor
import cromwell.core.{ExecutionEvent, CallOutputs}

/**
 * ADT representing the result of an execution of a BackendCall.
 */
sealed trait ExecutionResult

/**
 * A successful execution with resolved outputs.
 */
final case class SuccessfulExecution(outputs: CallOutputs,
                                     returnCode: Int,
                                     jobDetritusFiles: Map[String, Path],
                                     executionEvents: Seq[ExecutionEvent],
                                     resultsClonedFrom: Option[BackendJobDescriptor] = None) extends ExecutionResult

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
