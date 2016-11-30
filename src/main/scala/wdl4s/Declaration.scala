package wdl4s

import wdl4s.AstTools.EnhancedAstNode
import wdl4s.parser.WdlParser.{Ast, AstNode}
import wdl4s.types.{WdlArrayType, WdlType}

object DeclarationInterface {
  /**
    * Depending on who is asking, the type of a declaration can vary.
    * e.g
    * task a {
    *   ...
    *   output {
    *     String o
    *   }
    * }
    * workflow w {
    *   scatter (...) {
    *     call a
    *     String s = a.o # Inside the scatter a.o is a String
    *   }
    *
    *   Array[String] s2 = a.o # Outside the scatter it's an Array[String]
    * }
    */
  def relativeWdlType(from: Scope, target: DeclarationInterface, wdlType: WdlType): WdlType = {
    target.closestCommonAncestor(from) map { ancestor =>
      target.ancestrySafe.takeWhile(_ != ancestor).collect({
        case s: Scatter => s
      }).foldLeft(wdlType)((acc, _) => WdlArrayType(acc))
    } getOrElse wdlType
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
trait DeclarationInterface extends GraphNode {
  def wdlType: WdlType
  def expression: Option[WdlExpression]
  def ast: Ast

  def relativeWdlType(from: Scope): WdlType = DeclarationInterface.relativeWdlType(from, this, wdlType)

  def asTaskInput: Option[TaskInput] = expression match {
    case Some(expr) => None
    case None => Option(TaskInput(unqualifiedName, wdlType))
  }

  def asWorkflowInput: Option[WorkflowInput] = expression match {
    case Some(expr) => None
    case None => Some(WorkflowInput(fullyQualifiedName, wdlType))
  }

  def toWdlString: String = {
    val expr = expression.map(e => s" = ${e.toWdlString}").getOrElse("")
    s"${wdlType.toWdlString} $unqualifiedName$expr"
  }

  lazy val upstream: Set[GraphNode] = {
    val nodes = for {
      expr <- expression.toSeq
      variable <- expr.variableReferences
      node <- resolveVariable(variable.sourceString)
    } yield node
    nodes.toSet
  }

  lazy val downstream: Set[GraphNode] = {
    for {
      node <- namespace.descendants.collect({ 
        case n: GraphNode if n.fullyQualifiedName != fullyQualifiedName => n 
      })
      if node.upstream.contains(this)
    } yield node
  }

  override def toString: String = {
    s"[Declaration type=${wdlType.toWdlString} name=$unqualifiedName expr=${expression.map(_.toWdlString)}]"
  }
}

object Declaration {

  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter, parent: Option[Scope]): Declaration = {
    Declaration(
      ast.getAttribute("type").wdlType(wdlSyntaxErrorFormatter),
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
                       unqualifiedName: String,
                       expression: Option[WdlExpression],
                       override val parent: Option[Scope],
                       ast: Ast) extends DeclarationInterface
