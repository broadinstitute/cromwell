package wdl

import wdl.AstTools.EnhancedAstNode
import wdl4s.parser.WdlParser.Ast
import wom.types.WdlType

object WorkflowOutput {
  def apply(ast: Ast, syntaxErrorFormatter: WdlSyntaxErrorFormatter, parent: Option[Scope]): WorkflowOutput = {
    val wdlType = ast.getAttribute("type").wdlType(syntaxErrorFormatter)
    val name = ast.getAttribute("name").sourceString
    val expression = WdlExpression(ast.getAttribute("expression"))
    WorkflowOutput(name, wdlType, expression, ast, parent)
  }
}

case class WorkflowOutput(unqualifiedName: String, wdlType: WdlType, requiredExpression: WdlExpression, ast: Ast, override val parent: Option[Scope]) extends Output