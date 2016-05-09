package cromwell.engine.backend

import cromwell.core.CallOutputs
import cromwell.engine.ExecutionEventEntry

/**
 * ADT representing the result of an execution of a BackendCall.
 */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
sealed trait OldStyleExecutionResult
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case object OldStyleSuccessfulFinalCallExecution extends OldStyleExecutionResult

/**
 * A successful execution with resolved outputs.
 */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
final case class OldStyleSuccessfulBackendCallExecution(outputs: CallOutputs, executionEvents: Seq[ExecutionEventEntry], returnCode: Int, hash: ExecutionHash, resultsClonedFrom: Option[OldStyleBackendCallJobDescriptor] = None) extends OldStyleExecutionResult

/**
 * A user-requested abort of the command.
 */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case object OldStyleAbortedExecution extends OldStyleExecutionResult
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
sealed trait OldStyleFailedExecution extends OldStyleExecutionResult {
  def e: Throwable
  def returnCode: Option[Int]
}

/**
  * Failed execution, possibly having a return code.
  */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
final case class OldStyleNonRetryableFailedExecution(e: Throwable, returnCode: Option[Int] = None, events: Seq[ExecutionEventEntry] = Seq.empty) extends OldStyleFailedExecution
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
final case class OldStyleRetryableFailedExecution(e: Throwable, returnCode: Option[Int] = None, events: Seq[ExecutionEventEntry] = Seq.empty) extends OldStyleFailedExecution
