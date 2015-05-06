package cromwell.binding

import cromwell.binding.command.Command
import cromwell.binding.types.WdlType

/**
 * Represents a `task` declaration in a WDL file
 *
 * @param name Name of the task
 * @param command Abstract command defined in the `command` section
 * @param outputs Set of defined outputs in the `output` section of the task
 */
case class Task(name: String, command: Command, outputs: Set[TaskOutput]) {
  /**
   * Inputs to this task, as task-local names (i.e. not fully-qualified)
   *
   * @return Map of input name to type for that input
   */
  val inputs: Map[String, WdlType] = command.inputs

  override def toString: String = s"[Task name=$name command=$command]"
}
