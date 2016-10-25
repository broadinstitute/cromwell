package wdl4s

import wdl4s.AstTools.EnhancedAstNode
import wdl4s.types.WdlType
import wdl4s.parser.WdlParser.{Ast, AstNode}

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
trait DeclarationInterface extends Scope with GraphNode {
  def wdlType: WdlType
  def postfixQuantifier: Option[String]
  def expression: Option[WdlExpression]
  def ast: Ast

  def asTaskInput: Option[TaskInput] = expression match {
    case Some(expr) => None
    case None => Option(TaskInput(unqualifiedName, wdlType, postfixQuantifier))
  }

  def asWorkflowInput: Option[WorkflowInput] = expression match {
    case Some(expr) => None
    case None => Some(WorkflowInput(fullyQualifiedName, wdlType, postfixQuantifier))
  }

  def toWdlString: String = {
    val expr = expression.map(e => s" = ${e.toWdlString}").getOrElse("")
    s"${wdlType.toWdlString} $unqualifiedName$expr"
  }

  lazy val upstream: Set[Scope with GraphNode] = {
    val nodes = for {
      expr <- expression.toSeq
      variable <- expr.variableReferences
      node <- resolveVariable(variable.sourceString)
    } yield node
    nodes.toSet
  }

  lazy val downstream: Set[Scope with GraphNode] = {
    for {
      node <- namespace.descendants.collect({ 
        case n: GraphNode if n.fullyQualifiedName != fullyQualifiedName => n 
      })
      if node.upstream.contains(this)
    } yield node
  }

  override def toString(): String = {
    s"[Declaration name=$unqualifiedName expr=${expression.map(_.toWdlString)}]"
  }
}

object Declaration {
  val OptionalPostfixQuantifier = "?"
  
  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter, parent: Option[Scope]): Declaration = {
    Declaration(
      ast.getAttribute("type").wdlType(wdlSyntaxErrorFormatter),
      Option(ast.getAttribute("postfix")).map(_.sourceString),
      ast.getAttribute("name").sourceString,
      ast.getAttribute("expression") match {
        case a: AstNode => Option(WdlExpression(a))
        case _ => None
      },
      parent,
      ast
    )
  }
}

case class Declaration(wdlType: WdlType,
                       postfixQuantifier: Option[String],
                       unqualifiedName: String,
                       expression: Option[WdlExpression],
                       override val parent: Option[Scope],
                       ast: Ast) extends DeclarationInterface
