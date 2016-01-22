package wdl4s

import wdl4s.AstTools.EnhancedAstNode
import wdl4s.types.WdlType
import wdl4s.parser.WdlParser.{Ast, AstNode}

object Declaration {
  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Declaration = {
    UnscopedDeclaration(
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
trait Declaration {
  def wdlType: WdlType
  def postfixQuantifier: Option[String]
  def name: String
  def expression: Option[WdlExpression]
  def asTaskInput: Option[TaskInput] = expression match {
    case Some(expr) => None
    case None => Option(TaskInput(name, wdlType, postfixQuantifier))
  }

  def toWdlString: String = {
    val expr = expression.map(e => s" = ${e.toWdlString}").getOrElse("")
    s"${wdlType.toWdlString} $name$expr"
  }
}

case class UnscopedDeclaration(wdlType: WdlType, postfixQuantifier: Option[String], name: String, expression: Option[WdlExpression]) extends Declaration

object ScopedDeclaration {
  def apply(scope: Scope, decl: Declaration): ScopedDeclaration = {
    ScopedDeclaration(scope, decl.wdlType, decl.postfixQuantifier, decl.name, decl.expression)
  }
}

case class ScopedDeclaration(scope: Scope, wdlType: WdlType, postfixQuantifier: Option[String], name: String, expression: Option[WdlExpression]) extends Declaration {
  def fullyQualifiedName: FullyQualifiedName = s"${scope.fullyQualifiedName}.$name"
  def asWorkflowInput: Option[WorkflowInput] = expression match {
    case Some(expr) => None
    case None => Some(WorkflowInput(fullyQualifiedName, wdlType, postfixQuantifier))
  }
}