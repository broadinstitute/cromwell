package cromwell.binding

import cromwell.binding.types.WdlType
import cromwell.parser.WdlParser.{Terminal, Ast}

object Workflow {
  def apply(ast: Ast, calls: Seq[Call]): Workflow = {
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    new Workflow(name, calls)
  }
}

/**
 * Represents a `workflow` definition in WDL which currently
 * only supports a set of `call` declarations and a name for
 * the workflow
 *
 * @param name The name of the workflow
 * @param calls The set of `call` declarations
 */
case class Workflow(name: String, calls: Seq[Call]) extends Executable with Scope {
  calls foreach {c => c.setParent(this)}

  /** Parent node for this workflow.  Since we do not support nested
    * workflows currently, this is always `None`
    */
  val parent: Option[Scope] = None

  /**
   * All inputs for this workflow and their associated types.
   *
   * @return a Seq[WorkflowInput] representing the
   *         inputs that the user needs to provide to this workflow
   */
  def inputs: Seq[WorkflowInput] =
    for {call <- calls; input <- call.unsatisfiedInputs} yield input

  /**
   * All outputs for this workflow and their associated types
   *
   * @return a Map[FullyQualifiedName, WdlType] representing the union
   *         of all outputs from all `call`s within this workflow
   */
  def outputs: Map[FullyQualifiedName, WdlType] = {
    val outputs = for (call <- calls; output <- call.task.outputs) yield (s"${call.fullyQualifiedName}.${output.name}", output.wdlType)
    outputs.toMap
  }
}
