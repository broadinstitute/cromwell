package wdl4s

import wdl4s.AstTools.EnhancedAstNode
import wdl4s.types.WdlType
import wdl4s.parser.WdlParser.{Ast, AstNode}

object Declaration {
  def apply(ast: Ast, scopeFqn: FullyQualifiedName, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Declaration = {
    Declaration(
      scopeFqn,
      ast.getAttribute("type").wdlType(wdlSyntaxErrorFormatter),
      Option(ast.getAttribute("postfix")).map(_.sourceString),
      ast.getAttribute("name").sourceString,
      ast.getAttribute("expression") match {
        case a: AstNode => Some(WdlExpression(a))
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
case class Declaration(scopeFqn: FullyQualifiedName, wdlType: WdlType, postfixQuantifier: Option[String], name: String, expression: Option[WdlExpression]) {
  def asWorkflowInput: Option[WorkflowInput] = expression match {
    case Some(expr) => None
    case None => Some(WorkflowInput(fullyQualifiedName, wdlType, postfixQuantifier))
  }

  def asTaskInput: Option[TaskInput] = expression match {
    case Some(expr) => None
    case None => Some(TaskInput(name, wdlType, postfixQuantifier))
  }

  def fullyQualifiedName: FullyQualifiedName = s"$scopeFqn.$name"
}
