package wdl4s

import better.files.File
import wdl4s.AstTools.EnhancedAstNode
import wdl4s.parser.WdlParser.{Terminal, Ast, AstNode}

object Import {
  def apply(astNode: AstNode): Import = {
    val ast = astNode.asInstanceOf[Ast]
    val uri = ast.getAttribute("uri")
    val namespaceTerminal = Option(ast.getAttribute("namespace")).orElse(Option(uri)).map(_.asInstanceOf[Terminal]).get

    Import(uri.sourceString, namespaceTerminal)
  }
}

case class Import(uri: String, namespaceTerminal: Terminal) {
  val namespaceName: String = File(namespaceTerminal.sourceString).nameWithoutExtension
}
