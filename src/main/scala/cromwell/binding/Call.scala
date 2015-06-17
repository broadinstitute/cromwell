package cromwell.binding

import cromwell.binding.types.WdlType
import cromwell.binding.values.WdlValue
import cromwell.parser.AstTools

import scala.util.{Success, Try}

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
case class Call(alias: Option[String], taskFqn: FullyQualifiedName, task: Task, inputMappings: Map[String, WdlExpression], namespace: WdlNamespace) extends Scope {
  val name: String = alias getOrElse taskFqn

  private var _parent: Option[Scope] = None

  def parent: Option[Scope] = _parent

  def setParent(parent: Scope) = {
    if (this._parent.isEmpty) this._parent = Option(parent)
    else throw new UnsupportedOperationException("parent is write-once")
  }

  private def unsatisfiedTaskInputs: Seq[TaskInput] =
    task.inputs.filterNot {case i => inputMappings.contains(i.name)}

  /**
   * Map of task-local input names to type for every input that this
   * Call needs to execute the task
   *
   * @return Seq[WorkflowInput] representing each task-local input
   */
  def unsatisfiedInputs: Seq[WorkflowInput] = {
    unsatisfiedTaskInputs.map {i =>
      WorkflowInput(s"$fullyQualifiedName.${i.name}", i.types, i.postfixQuantifier)
    }
  }

  override def toString: String = s"[Call name=$name, task=$task]"

  /**
   * Find all calls upon which this call immediately depends, i.e. the result of this
   * does not include transitive dependencies.  Currently this only works for member
   * access expressions with a literal LHS, e.g.:
   *
   * {{{
   *   call cgrep {
   *     input: in_file=ps.procs
   *   }
   * }}}
   *
   * Here `ps` would be the prerequisite call for `cgrep`.
   *
   * Calls are de-duplicated into a returned `Set`, so if one call expresses dependencies on
   * a prerequisite call multiple times (i.e. has multiple inputs depending on multiple outputs
   * of a prerequisite call), that prerequisite call will appear only once in the output.
   */
  def prerequisiteCalls(): Iterable[Call] = {
    for {
      expr <- inputMappings.values
      ast <- AstTools.findTopLevelMemberAccesses(expr.ast)
      call <- namespace.getCallFromMemberAccessAst(ast) match {
        case Success(c:Call) => Vector(c)
        case _ => Vector()
      }
    } yield call
  }

  /**
   * Instantiate the abstract command line corresponding to this call using the specified inputs.
   */
  def instantiateCommandLine(inputs: CallInputs): Try[String] = task.command.instantiate(inputs)

  /**
   * Return the docker configuration value associated with this `Call`, if any.
   */
  def docker: Option[String] = task.runtimeAttributes.docker
}
