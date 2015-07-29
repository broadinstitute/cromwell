package cromwell.binding

import cromwell.binding.types.WdlType
import AstTools.EnhancedAstNode
import cromwell.parser.WdlParser.Ast

object TaskOutput {
  def apply(ast: Ast, syntaxErrorFormatter: WdlSyntaxErrorFormatter): TaskOutput = {
    val wdlType = ast.getAttribute("type").wdlType(syntaxErrorFormatter)
    val name = ast.getAttribute("var").sourceString()
    val expression = ast.getAttribute("expression")
    new TaskOutput(name, wdlType, new WdlExpression(expression))
  }
}

case class TaskOutput(name: String, wdlType: WdlType, expression: WdlExpression)
