package cromwell.engine.backend

import cromwell.binding.WdlExpression.ScopedLookupFunction
import cromwell.binding._
import cromwell.binding.values.WdlValue

import scala.util.Try

/**
 * Trait to be implemented by concrete backends.
 */
trait Backend {

  /**
   * Return a possibly altered copy of inputs reflecting any localization of input file paths that might have
   * been performed for this `Backend` implementation.
   */
  def adjustInputPaths(call: Call, inputs: CallInputs): CallInputs

  /**
   * Do whatever work is required to initialize the workflow, returning a copy of
   * the coerced inputs present in the `WorkflowDescriptor` with any input `WdlFile`s
   * adjusted for the host workflow execution path.
   */
  def initializeForWorkflow(workflow: WorkflowDescriptor): HostInputs

  /**
   * Execute the specified command line using the provided symbol store, evaluating the task outputs to produce
   * a mapping of local task output names to WDL values.
   */
  def executeCommand(commandLine: String, workflowDescriptor: WorkflowDescriptor, call: Call, scopedLookupFunction: ScopedLookupFunction): Map[String, Try[WdlValue]]

}
