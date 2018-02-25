package wdl.draft2.model

import better.files.File
import wdl.draft2.model.AstTools.EnhancedAstNode
import wdl.draft2.parser.WdlParser.{Ast, AstNode, Terminal}

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
