package cromwell.engine.backend

import cromwell.engine.{CallEngineFunctions, CallOutputs, ExecutionEventEntry, ExecutionHash}
import wdl4s._
import wdl4s.values._

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

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

trait BackendCall {
  /**
   * The Workflow and Call to invoke.  It is assumed that in the creation
   * of a BackendCall object that the 'call' would be within the workflow
   */
  def jobDescriptor: BackendCallJobDescriptor
  def workflowDescriptor = jobDescriptor.workflowDescriptor
  def key = jobDescriptor.key
  def call = key.scope

  /**
    * Inputs to the call.  For example, if a call's task specifies a command like this:
    *
    * File some_dir
    * command {
    *   ls ${some_dir}
    * }
    *
    * Then locallyQualifiedInputs could be Map("some_dir" -> WdlFile("/some/path"))
    */
  def locallyQualifiedInputs: CallInputs = jobDescriptor.locallyQualifiedInputs

  /**
   * Backend which will be used to execute the Call
   */
  def backend: Backend

  def callRootPathWithBaseRoot(baseRoot: String) = jobDescriptor.callRootPathWithBaseRoot(baseRoot)

  def callRootPath = jobDescriptor.callRootPath

  /**
    * Attempt to evaluate all the ${...} tags in a command and return a String representation
    * of the command.  This could fail for a variety of reasons related to expression evaluation
    * which is why it returns a Try[String]
    */
  def instantiateCommand: Try[String] = jobDescriptor.instantiateCommand

  /**
   * Implementation of CallEngineFunctions, used to evaluate WdlExpressions
   */
  def callEngineFunctions: CallEngineFunctions = jobDescriptor.callEngineFunctions

  /**
   * Function used to get the value for identifiers in expressions.  For example, the
   * expression `read_lines(my_file_var)` would have to call lookupFunction()("my_file_var")
   * during expression evaluation
   */
  def lookupFunction(evaluatedValues: Map[String, WdlValue]): String => WdlValue = jobDescriptor.lookupFunction(evaluatedValues)

  /** Initiate execution, callers can invoke `poll` once this `Future` completes successfully. */
  def execute(implicit ec: ExecutionContext): Future[ExecutionHandle]

  /**
   * The default implementation of this method is not expected to be called and simply throws an `NotImplementedError`.
   * If the corresponding backend does not override `Backend#findResumableExecutions` to return resumable executions,
   * this method will not be called.  If the backend does override `Backend#findResumableExecutions`, the corresponding
   * `BackendCall` should override this method to actually do the resumption work.
   */
  def resume(jobKey: JobKey)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    throw new NotImplementedError(s"resume() called on a non-resumable BackendCall: $this")
  }

  def useCachedCall(cachedBackendCall: BackendCall)(implicit ec: ExecutionContext): Future[ExecutionHandle] = jobDescriptor.useCachedCall(cachedBackendCall.jobDescriptor)

  /**
    * Return CallLogs which contains the stdout/stderr of the particular call
    */
  def stdoutStderr: CallLogs

  @throws[IllegalArgumentException]
  lazy val runtimeAttributes = jobDescriptor.callRuntimeAttributes

  /**
    * Compute a hash that uniquely identifies this call
    */
  def hash(implicit ec: ExecutionContext): Future[ExecutionHash] = jobDescriptor.hash

  /**
   * Using the execution handle from the previous execution, resumption, or polling attempt, poll the execution
   * of this `BackendCall`.
   */
  def poll(previous: ExecutionHandle)(implicit ec: ExecutionContext): Future[ExecutionHandle] = jobDescriptor.poll(previous)

}
