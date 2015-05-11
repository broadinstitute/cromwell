package cromwell.binding

import cromwell.binding.types.WdlType

/**
 * Represents a `call` block in a WDL workflow.  Calls wrap tasks
 * and optionally provide a subset of the inputs for the task (as inputMappings member).
 * All remaining inputs to the task that are not declared in the `input` section
 * of the `call` block are called unsatisfiedInputs
 *
 * @param alias The alias for this call.  If two calls in a workflow use the same task
 *              then one of them needs to be aliased to a different name
 * @param task The task that this `call` will invoke
 * @param inputMappings A map of task-local input names and corresponding expression for the
 *                      value of those inputs
 * @todo Validate that the keys in inputMappings correspond to actual parameters in the task
 */
case class Call(alias: Option[String], task: Task, inputMappings: Map[String, WdlExpression]) extends Scope {
  val name: String = alias getOrElse task.name

  private var _parent: Option[Scope] = None

  def parent: Option[Scope] = _parent

  def setParent(parent: Scope) = {
    if (this._parent == None) this._parent = Option(parent)
    else throw new UnsupportedOperationException("parent is write-once")
  }

  /**
   * Map of task-local input names to type for every input that this
   * Call needs to execute the task
   *
   * @return Map[String, WdlType] representing each task-local input
   */
  def unsatisfiedInputs: Map[String, WdlType] = task.inputs.filterNot { case (k, v) => inputMappings.contains(k) }
  
  override def toString: String = s"[Call name=$name, task=$task]"
}
