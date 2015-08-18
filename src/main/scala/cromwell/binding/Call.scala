package cromwell.binding

import cromwell.binding.AstTools.EnhancedAstNode
import cromwell.parser.WdlParser.{Ast, SyntaxError, Terminal}

import scala.util.{Success, Try}
import scala.language.postfixOps

object Call {
  def apply(ast: Ast,
            namespaces: Seq[WdlNamespace],
            tasks: Seq[Task],
            wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Call = {
    val alias: Option[String] = ast.getAttribute("alias") match {
      case x: Terminal => Option(x.getSourceString)
      case _ => None
    }

    val taskName = ast.getAttribute("task").sourceString
    val task = WdlNamespace.findTask(taskName, namespaces, tasks) getOrElse {
      throw new SyntaxError(wdlSyntaxErrorFormatter.callReferencesBadTaskName(ast, taskName))
    }

    val callInputSectionMappings = processCallInput(ast, wdlSyntaxErrorFormatter)

    callInputSectionMappings foreach { case (taskParamName, expression) =>
      task.inputs find { taskInput => taskInput.name == taskParamName } getOrElse {
        /*
          FIXME-ish
          It took me a while to figure out why this next part is necessary and it's kind of hokey.
          All the syntax error formatting requires ASTs and this is a way to get the input's AST back.
          Perhaps an intermediary object that contains the AST as well and then the normal Map in the Call?

          FIXME: It'd be good to break this whole thing into smaller chunks
         */
        val callInput = AstTools.callInputSectionIOMappings(ast, wdlSyntaxErrorFormatter) find {
          _.getAttribute("key").sourceString == taskParamName
        } getOrElse {
          throw new SyntaxError(s"Can't find call input: $taskParamName")
        }
        throw new SyntaxError(wdlSyntaxErrorFormatter.callReferencesBadTaskInput(callInput, task.ast))
      }
    }

    new Call(alias, taskName, task, callInputSectionMappings)
  }

  private def processCallInput(ast: Ast,
                               wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Map[String, WdlExpression] = {
    AstTools.callInputSectionIOMappings(ast, wdlSyntaxErrorFormatter) map { a =>
      val key = a.getAttribute("key").sourceString
      val expression = new WdlExpression(a.getAttribute("value"))
      (key, expression)
    } toMap
  }
}

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
 */
case class Call(alias: Option[String],
                taskFqn: FullyQualifiedName,
                task: Task,
                inputMappings: Map[String, WdlExpression]) extends Scope {
  val name: String = alias getOrElse taskFqn

  /*
    TODO/FIXME: Since a Workflow's Calls *must* have a parent a better way to handle this would be to use an ADT
    where one type has a parent and one does not, and the Workflow can only take the former. I went down that
    road a bit but Scope requires parent to be an Option[Scope] so it's turtles all the way down (well, up in this case)
   */
  private var _parent: Option[Scope] = None

  def parent: Option[Scope] = _parent

  def setParent(parent: Scope) = {
    if (this._parent.isEmpty) this._parent = Option(parent)
    else throw new UnsupportedOperationException("parent is write-once")
  }

  private def unsatisfiedTaskInputs: Seq[TaskInput] = task.inputs.filterNot {case i => inputMappings.contains(i.name)}

  /**
   * Returns a Seq[WorkflowInput] representing the inputs to the call that are
   * needed before its command can be constructed. This excludes inputs that
   * are satisfied via the 'input' section of the Call definition.
   */
  def unsatisfiedInputs: Seq[WorkflowInput] = unsatisfiedTaskInputs.map {i =>
    WorkflowInput(s"$fullyQualifiedName.${i.name}", i.wdlType, i.postfixQuantifier)
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
  /*
    TODO/FIXME: Not happy w/ having to include the namespace but since Calls are now built before there's a namespace the Call loses it's WdlNamespace field.
   */
  def prerequisiteCalls(namespace: NamespaceWithWorkflow): Iterable[Call] = {
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

  /**
   * Whether the call should be considered to have failed if any stderr is generated.
   */
  def failOnStderr: Boolean = task.runtimeAttributes.failOnStderr
}
