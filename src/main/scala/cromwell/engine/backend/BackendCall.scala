package cromwell.engine.backend

import cromwell.binding._
import cromwell.binding.expression.WdlStandardLibraryFunctions
import cromwell.binding.values.WdlValue
import cromwell.engine.workflow.CallKey

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

trait BackendCall {

  /**
   * The Workflow and Call to invoke.  It is assumed that in the creation
   * of a BackendCall object that the 'call' would be within the workflow
   */
  def workflowDescriptor: WorkflowDescriptor
  def key: CallKey
  def call = key.scope

  /**
   * Backend which will be used to execute the Call
   */
  def backend: Backend

  /**
   * Inputs to the call.  For example, if a call's task specifies a command like this:
   *
   * command {
   *   File some_dir
   *   ls ${some_dir}
   * }
   *
   * Then locallyQualifiedInputs could be Map("some_dir" -> WdlFile("/some/path"))
   */
  def locallyQualifiedInputs: CallInputs

  /**
   * Implementation of the WDL Standard Library Functions, used to evaluate WdlExpressions
   */
  def engineFunctions: WdlStandardLibraryFunctions

  /**
   * Function used to get the value for identifiers in expressions.  For example, the
   * expression `read_lines(my_file_var)` would have to call lookupFunction()("my_file_var")
   * during expression evaluation
   */
  def lookupFunction: String => WdlValue

  /**
   * Attempt to evaluate all the ${...} tags in a command and return a String representation
   * of the command.  This could fail for a variety of reasons related to expression evaluation
   * which is why it returns a Try[String]
   */
  def instantiateCommand: Try[String] = {
    val backendInputs = backend.adjustInputPaths(call, locallyQualifiedInputs)
    call.instantiateCommandLine(backendInputs, engineFunctions)
  }

  /**
   * Block while executing this call.  Return a Success(Map[CallOutputs]) if the
   * Call ran successfully.  Otherwise return a Failure.
   */
  def execute: ExecutionResult
}
