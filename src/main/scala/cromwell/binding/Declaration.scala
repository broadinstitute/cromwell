package cromwell.binding

import cromwell.binding.AstTools.EnhancedAstNode
import cromwell.binding.types.WdlType
import cromwell.parser.WdlParser.Ast

object Declaration {
  def apply(ast: Ast, scopeFqn: FullyQualifiedName, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Declaration = {
    Declaration(
      scopeFqn,
      ast.getAttribute("type").wdlType(wdlSyntaxErrorFormatter),
      ast.getAttribute("name").sourceString,
      ast.getAttribute("expression") match {
        case a: Ast => Some(WdlExpression(a))
        case _ => None
      }
    )
  }
}

/**
 * Represents a declaration which can show up in a workflow or a task context.  For example
 *
 * task test {
 *   File test_file
 *   command { ... }
 * }
 *
 * workflow wf {
 *   String wf_string = "abc"
 *   call test { input: s=wf_string }
 * }
 *
 * Both the definition of test_file and wf_string are declarations
 */
case class Declaration(scopeFqn: FullyQualifiedName, wdlType: WdlType, name: String, expression: Option[WdlExpression]) {
  def asWorkflowInput: Option[WorkflowInput] = expression match {
    case Some(expr) => None
    case None => Some(WorkflowInput(fullyQualifiedName, wdlType, postfixQuantifier = None))
  }

  def asTaskInput: Option[TaskInput] = expression match {
    case Some(expr) => None
    case None => Some(TaskInput(name, wdlType, postfixQuantifier = None))
  }

  def fullyQualifiedName: FullyQualifiedName = s"$scopeFqn.$name"
}
