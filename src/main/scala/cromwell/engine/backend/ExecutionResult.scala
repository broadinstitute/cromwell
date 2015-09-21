package cromwell.engine.backend

import cromwell.binding._

import scala.util.Try


/**
 * ADT representing the result of an execution of a BackendCall.
 */
sealed trait ExecutionResult

/**
 * A successful execution with resolved outputs.
 */
final case class SuccessfulExecution(outputs: CallOutputs) extends ExecutionResult

/**
 * A user-requested abort of the command.
 */
case object AbortedExecution extends ExecutionResult

/**
 * Failed execution, possibly having a return code.
 */
final case class FailedExecution(e: Throwable, returnCode: Option[Int] = None) extends ExecutionResult

