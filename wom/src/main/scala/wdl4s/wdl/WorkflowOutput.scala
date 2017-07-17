package wdl4s.wdl

import wdl4s.wdl.AstTools.EnhancedAstNode
import wdl4s.parser.WdlParser.Ast
import wdl4s.wdl.types.WdlType

object WorkflowOutput {
  def apply(ast: Ast, syntaxErrorFormatter: WdlSyntaxErrorFormatter, parent: Option[Scope]): WorkflowOutput = {
    val wdlType = ast.getAttribute("type").wdlType(syntaxErrorFormatter)
    val name = ast.getAttribute("name").sourceString
    val expression = WdlExpression(ast.getAttribute("expression"))
    WorkflowOutput(name, wdlType, expression, ast, parent)
  }
}

case class WorkflowOutput(unqualifiedName: String, wdlType: WdlType, requiredExpression: WdlExpression, ast: Ast, override val parent: Option[Scope]) extends Output