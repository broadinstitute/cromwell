package cromwell.binding

import cromwell.binding.AstTools.EnhancedAstNode
import cromwell.binding.expression.WdlFunctions
import cromwell.binding.values.WdlValue
import cromwell.parser.WdlParser.{Ast, SyntaxError, Terminal}
import scala.util.Try
import scala.language.postfixOps

object Call {
  def apply(ast: Ast,
            namespaces: Seq[WdlNamespace],
            tasks: Seq[Task],
            wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter,
            parent: Option[Scope]): Call = {
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

    val prerequisiteCallNames = callInputSectionMappings.values flatMap { _.prerequisiteCallNames } toSet

    new Call(alias, taskName, task, prerequisiteCallNames, callInputSectionMappings, parent)
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
                prerequisiteCallNames: Set[LocallyQualifiedName],
                inputMappings: Map[String, WdlExpression],
                parent: Option[Scope]) extends Scope {
  val unqualifiedName: String = alias getOrElse taskFqn

  override lazy val prerequisiteScopes: Set[Scope] = {
    val parent = this.parent.get // FIXME: In a world where Call knows it has a parent this wouldn't be icky
    val parentPrereq = if (parent == this.rootWorkflow) Nil else Set(parent)
    prerequisiteCalls ++ parentPrereq
  }

  /**
   * Returns a Seq[WorkflowInput] representing the inputs to the call that are
   * needed before its command can be constructed. This excludes inputs that
   * are satisfied via the 'input' section of the Call definition.
   */
  def unsatisfiedInputs: Seq[WorkflowInput] = for {
    i <- task.inputs if !inputMappings.contains(i.name)
  } yield WorkflowInput(s"$fullyQualifiedName.${i.name}", i.wdlType, i.postfixQuantifier)

  override def toString: String = s"[Call name=$unqualifiedName, task=$task]"

  /**
   * Instantiate the abstract command line corresponding to this call using the specified inputs.
   */
  def instantiateCommandLine(inputs: CallInputs, functions: WdlFunctions[WdlValue]): Try[String] =
    task.instantiateCommand(inputs, functions)

  /**
   * Return the docker configuration value associated with this `Call`, if any.
   */
  def docker: Option[String] = task.runtimeAttributes.docker

  /**
   * Whether the call should be considered to have failed if any stderr is generated.
   */
  def failOnStderr: Boolean = task.runtimeAttributes.failOnStderr

  /**
   * Whether the call should not be considered to have failed if a non-zero result code is generated,
   * or one of the supplied result code is generated.
   */
  def continueOnReturnCode: ContinueOnReturnCode = task.runtimeAttributes.continueOnReturnCode

  override def rootWorkflow: Workflow = parent map { _.rootWorkflow } getOrElse { throw new IllegalStateException("Call not in workflow") }
}
