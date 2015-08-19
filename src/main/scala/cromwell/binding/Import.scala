package cromwell.binding

import cromwell.binding.AstTools.EnhancedAstNode
import cromwell.parser.WdlParser.{Ast, AstNode}

object Import {
  def apply(astNode: AstNode): Import = {
    val uri = astNode.asInstanceOf[Ast].getAttribute("uri").sourceString
    val importNamespace = Option(astNode.asInstanceOf[Ast].getAttribute("namespace"))
    Import(uri, importNamespace)
  }
}

// FIXME: I dislike dragging the AST along but it's necessary for "compile" time error syntax highlighting, argh
case class Import(uri: String, namespaceAst: Option[AstNode]) {
  val namespace: Option[String] = namespaceAst map { _.sourceString }
}
