package cromwell.engine.backend

import cromwell.binding.WdlStandardLibraryFunctions

/**
 * TaskExecutionContext represents the static information about a specific Call's invocation.
 *
 * For the LocalBackend, this would hold information like the Call's current working directory
 * where it executes the command and a path to stdout/stderr files (which might not exist yet).
 *
 * An object with this trait is available to all WDL standard library functions.
 * This makes it possible to implement WDL functions like read_lines and write_lines,
 * which need have some execution context in order to work.
 */
trait TaskExecutionContext {

  /**
   * Return the implementation of all WDL functions (e.g. read_lines(), stdout()) used
   * by this backend.  One situation in which this is called is to get the function implementations
   * needed to instantiate the command line, which might contain a function call in it.
   *
   * See setupCallEnvironment() or TaskExecutionContext for more details
   *
   * For example: ./script ${write_lines(some_array)}
   */
  def engineFunctions: WdlStandardLibraryFunctions
}
